import json
import logging
import os
import re
import shutil
from collections import defaultdict
from datetime import date, datetime, timezone
from multiprocessing.pool import ThreadPool
from pathlib import Path
from pickle import dump
from pickle import load
from urllib.parse import unquote, urlparse
from urllib.request import Request, urlopen

import subprocess
from Bio import SeqIO
from tqdm import tqdm  # pip3 install tqdm

logger = logging.getLogger('update_mlst')
REST_BASE = 'https://rest.pubmlst.org'
DBASES_COMPAT_FILE = 'dbases.xml'
SCHEME_LIST_FILE = 'scheme_list.pkl'
SCHEME_CATALOG_FILE = 'scheme_catalog.json'
SCHEME_CATALOG_META_FILE = 'scheme_catalog_meta.json'
SCHEME_LOCAL_META_FILE = '.fastmlst_scheme_meta.json'
# BIGSdb sources whose allele FASTA must not be merged into mlst.fasta (BLAST index).
# rMLST is ribosomal and must not be mixed with standard MLST allele panels.
SCHEME_DATABASES_EXCLUDED_FROM_MLST_FASTA = ('pubmlst_rmlst_seqdef',)


def scheme_excluded_from_mlst_concat(codename):
    """True if alleles for this scheme folder must be omitted from the shared mlst.fasta."""
    cn = (codename or '').strip()
    if not cn:
        return False
    for db in SCHEME_DATABASES_EXCLUDED_FROM_MLST_FASTA:
        if cn == db or cn.startswith(f'{db}_'):
            return True
    return False


def _package_bundle_dir():
    """Directory shipped with the fastmlst package (read-only fallback catalog)."""
    return Path(__file__).resolve().parent / 'bundle'
AUTH_TEST_DB = 'pubmlst_neisseria_seqdef'
PUBMLST_AUTHORIZE_URL = 'https://pubmlst.org/bigsdb'
CONFIG_DIR = Path.home() / '.config' / 'fastmlst'
CONFIG_FILE = CONFIG_DIR / 'pubmlst_oauth.json'
TOKENS_FILE = CONFIG_DIR / 'pubmlst_tokens.json'

# NEW: Define a function to override the default database location
def set_pathdb(custom_path):
    """
    Override the default PubMLST database directory.

    Parameters:
        custom_path (str): Custom path for the MLST database directory.
    """
    global pathdb
    pathdb = Path(custom_path)
    pathdb.mkdir(parents=True, exist_ok=True)

# Define the path to the .cache directory in the user's home directory
home_dir = Path.home()
cache_dir = home_dir / '.cache' / 'fastmlst'

# Ensure the cache directory exists
cache_dir.mkdir(parents=True, exist_ok=True)

# Default pathdb points to the pubmlst folder in the default cache
pathdb = cache_dir / 'pubmlst'

necessary_file = ['mlst.fasta.nhr', 'mlst.fasta.nsq', 'mlst.fasta.nin']
_oauth_credentials = None


def save_obj(obj, name):
    with open(name, 'wb') as f:
        dump(obj, f, 2)


def load_obj(name):
    with open(name, 'rb') as f:
        return load(f)


def _ensure_config_dir():
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    os.chmod(str(CONFIG_DIR), 0o700)


def save_pubmlst_client_credentials(client_id, client_secret):
    if not client_id or not client_secret:
        raise ValueError('client_id and client_secret are required.')
    _ensure_config_dir()
    payload = {'client_id': client_id, 'client_secret': client_secret}
    CONFIG_FILE.write_text(json.dumps(payload, indent=2), encoding='utf-8')
    os.chmod(str(CONFIG_FILE), 0o600)


def load_pubmlst_client_credentials():
    if not CONFIG_FILE.is_file():
        return None, None
    try:
        payload = json.loads(CONFIG_FILE.read_text(encoding='utf-8'))
    except Exception:
        return None, None
    return payload.get('client_id'), payload.get('client_secret')


def set_oauth_credentials(client_id=None, client_secret=None,
                          access_token=None, access_secret=None,
                          verifier=None):
    global _oauth_credentials
    if not client_id or not client_secret:
        stored_id, stored_secret = load_pubmlst_client_credentials()
        client_id = client_id or stored_id
        client_secret = client_secret or stored_secret
    values = [client_id, client_secret, access_token, access_secret, verifier]
    if any(values):
        if not client_id or not client_secret:
            raise ValueError('OAuth configuration requires client_id and client_secret.')
        if bool(access_token) != bool(access_secret):
            raise ValueError('If access_token is provided, access_secret is also required.')
        _oauth_credentials = {'client_id': client_id, 'client_secret': client_secret}
        if access_token and access_secret:
            _oauth_credentials['access_token'] = access_token
            _oauth_credentials['access_secret'] = access_secret
        if verifier:
            _oauth_credentials['verifier'] = verifier
    else:
        _oauth_credentials = None


def _oauth_cache_path():
    _ensure_config_dir()
    return TOKENS_FILE


def _load_cached_access_token(client_id):
    cache = _oauth_cache_path()
    if not cache.is_file():
        return None, None
    try:
        payload = json.loads(cache.read_text(encoding='utf-8'))
    except Exception:
        return None, None
    record = payload.get(client_id, {})
    return record.get('access_token'), record.get('access_secret')


def _save_cached_access_token(client_id, access_token, access_secret):
    cache = _oauth_cache_path()
    payload = {}
    if cache.is_file():
        try:
            payload = json.loads(cache.read_text(encoding='utf-8'))
        except Exception:
            payload = {}
    payload[client_id] = {'access_token': access_token, 'access_secret': access_secret}
    cache.write_text(json.dumps(payload, indent=2), encoding='utf-8')
    os.chmod(str(cache), 0o600)


def connect_pubmlst(client_id=None, client_secret=None, verifier=None):
    set_oauth_credentials(client_id=client_id, client_secret=client_secret, verifier=verifier)
    if not _oauth_credentials:
        raise RuntimeError('PubMLST credentials are missing. Provide client_id/client_secret.')
    save_pubmlst_client_credentials(_oauth_credentials['client_id'], _oauth_credentials['client_secret'])
    _build_authenticated_session()


def _exchange_client_credentials_for_access_tokens(client_id, client_secret, verifier=None):
    import requests
    from requests_oauthlib import OAuth1

    request_token_url = f'{REST_BASE}/db/{AUTH_TEST_DB}/oauth/get_request_token'
    access_token_url = f'{REST_BASE}/db/{AUTH_TEST_DB}/oauth/get_access_token'

    req_response = requests.get(
        request_token_url,
        auth=OAuth1(
            client_id,
            client_secret=client_secret,
            callback_uri='oob',
            signature_type='query'
        ),
        headers={'User-Agent': 'FastMLST'}
    )
    if req_response.status_code != 200:
        raise RuntimeError(
            f'Unable to obtain OAuth request token from PubMLST (HTTP {req_response.status_code}): '
            f'{req_response.text}'
        )
    req_tokens = req_response.json()
    request_token = req_tokens.get('oauth_token')
    request_secret = req_tokens.get('oauth_token_secret')
    if not request_token or not request_secret:
        raise RuntimeError('Unable to obtain OAuth request token from PubMLST.')

    auth_url = f'{PUBMLST_AUTHORIZE_URL}?db={AUTH_TEST_DB}&page=authorizeClient&oauth_token={request_token}'
    if not verifier:
        print('\nAuthorize FastMLST in your browser and paste the verifier code:')
        print(auth_url)
        verifier = input('OAuth verifier: ').strip()
    if not verifier:
        raise RuntimeError('OAuth verifier was not provided.')

    access_response = requests.get(
        access_token_url,
        auth=OAuth1(
            client_id,
            client_secret=client_secret,
            resource_owner_key=request_token,
            resource_owner_secret=request_secret,
            verifier=verifier,
            signature_type='query'
        ),
        headers={'User-Agent': 'FastMLST'}
    )
    if access_response.status_code != 200:
        raise RuntimeError(
            f'Unable to obtain OAuth access token from PubMLST (HTTP {access_response.status_code}): '
            f'{access_response.text}'
        )
    access_tokens = access_response.json()
    access_token = access_tokens.get('oauth_token')
    access_secret = access_tokens.get('oauth_token_secret')
    if not access_token or not access_secret:
        raise RuntimeError('Unable to obtain OAuth access token from PubMLST.')
    return access_token, access_secret


def _build_authenticated_session():
    if not _oauth_credentials:
        return None
    try:
        import requests
        from requests_oauthlib import OAuth1Session
    except ImportError as exc:
        raise RuntimeError(
            'OAuth credentials were provided but required dependencies are missing. '
            'Install requests and requests-oauthlib.'
        ) from exc

    access_token = _oauth_credentials.get('access_token')
    access_secret = _oauth_credentials.get('access_secret')
    if not access_token or not access_secret:
        cached_token, cached_secret = _load_cached_access_token(_oauth_credentials['client_id'])
        if cached_token and cached_secret:
            access_token, access_secret = cached_token, cached_secret
        else:
            access_token, access_secret = _exchange_client_credentials_for_access_tokens(
                _oauth_credentials['client_id'],
                _oauth_credentials['client_secret'],
                verifier=_oauth_credentials.get('verifier')
            )
            _save_cached_access_token(_oauth_credentials['client_id'], access_token, access_secret)

    session_bootstrap = OAuth1Session(
        client_key=_oauth_credentials['client_id'],
        client_secret=_oauth_credentials['client_secret'],
        resource_owner_key=access_token,
        resource_owner_secret=access_secret,
        signature_type='query',
    )
    token_url = f'{REST_BASE}/db/{AUTH_TEST_DB}/oauth/get_session_token'
    token_response = session_bootstrap.get(token_url, headers={'User-Agent': 'FastMLST'})
    if token_response.status_code != 200:
        raise RuntimeError(
            f'Unable to obtain PubMLST session token (HTTP {token_response.status_code}). '
            'Verify your OAuth credentials and that your account is registered.'
        )
    payload = token_response.json()
    session_token = payload.get('oauth_token')
    session_secret = payload.get('oauth_token_secret')
    if not session_token or not session_secret:
        raise RuntimeError('Invalid PubMLST OAuth session token response.')

    session = OAuth1Session(
        client_key=_oauth_credentials['client_id'],
        client_secret=_oauth_credentials['client_secret'],
        resource_owner_key=session_token,
        resource_owner_secret=session_secret,
        signature_type='query',
    )
    session.headers.update({'User-Agent': 'FastMLST'})
    # keep reference for isinstance checks without optional imports elsewhere
    session._fastmlst_requests_module = requests
    return session


def fetch_json(url, session=None):
    if session is not None:
        response = session.get(url, headers={'Accept': 'application/json'})
        response.raise_for_status()
        return response.json()
    request = Request(url, headers={'Accept': 'application/json', 'User-Agent': 'FastMLST'})
    with urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode('utf-8'))


def fetch_text(url, content_type='text/plain', session=None):
    if session is not None:
        response = session.get(url, headers={'Accept': content_type})
        response.raise_for_status()
        return response.text
    request = Request(url, headers={'Accept': content_type, 'User-Agent': 'FastMLST'})
    with urlopen(request, timeout=120) as response:
        return response.read().decode('utf-8')


def _text_nbytes(text):
    if not text:
        return 0
    return len(text.encode('utf-8'))


def scheme_storage_codename(database, scheme_id):
    """
    Directory name under schemes/: one folder per PubMLST database + scheme id.
    Avoids collisions from generic API descriptions (e.g. 'Extended MLST' -> 'extended').
    """
    db = (database or '').strip() or 'unknown_db'
    sid = int(scheme_id)
    safe = ''.join(c if (c.isalnum() or c == '_') else '_' for c in db)
    safe = safe.strip('_') or 'scheme'
    return f'{safe}_{sid}'


def ensure_scheme_codename(database, scheme_id, description=''):
    """
    Stable unique codename for REST-sourced schemes. ``description`` is ignored;
    kept for call-site compatibility. Never empty (see scheme_storage_codename).
    """
    return scheme_storage_codename(database, scheme_id)


def infer_locus_name_from_url(locus_url):
    parsed = urlparse(locus_url)
    return unquote(parsed.path.rstrip('/').split('/')[-1])


def remote_version_from_scheme_detail(detail):
    """Subset of GET .../schemes/{id} used to detect remote changes."""
    if not detail:
        return {}
    return {
        'scheme_id': detail.get('id'),
        'last_updated': detail.get('last_updated'),
        'last_added': detail.get('last_added'),
        'records': detail.get('records'),
        'locus_count': detail.get('locus_count'),
    }


def fetch_scheme_remote_version(database, scheme_id, session=None):
    scheme_url = f'{REST_BASE}/db/{database}/schemes/{scheme_id}'
    detail = fetch_json(scheme_url, session=session)
    return remote_version_from_scheme_detail(detail)


def local_scheme_meta_path(codename):
    return pathdb / 'schemes' / codename / SCHEME_LOCAL_META_FILE


def load_local_scheme_meta(codename):
    p = local_scheme_meta_path(codename)
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text(encoding='utf-8'))
    except Exception:
        return None


def save_local_scheme_meta(codename, database, scheme_id, description, remote_snapshot):
    p = local_scheme_meta_path(codename)
    p.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        'database': database,
        'scheme_id': scheme_id,
        'description': description,
        'remote_snapshot': remote_snapshot,
        'saved_at': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
    }
    p.write_text(json.dumps(payload, indent=2), encoding='utf-8')


def local_scheme_locus_stems(codename):
    d = pathdb / 'schemes' / codename
    if not d.is_dir():
        return []
    return sorted(p.stem for p in d.glob('*.tfa'))


def local_scheme_matches_remote(codename, database, scheme_id, session=None):
    """
    True if local meta and allele files exist and remote API reports the same
    version fingerprint (records / locus_count / last_updated when available).
    """
    stems = local_scheme_locus_stems(codename)
    if not stems:
        return False
    profiles = pathdb / 'schemes' / codename / f'{codename}.txt'
    if not profiles.is_file():
        return False
    meta = load_local_scheme_meta(codename)
    if not meta:
        return False
    try:
        meta_sid = int(meta.get('scheme_id'))
    except (TypeError, ValueError):
        return False
    if meta.get('database') != database or meta_sid != int(scheme_id):
        return False
    try:
        remote = fetch_scheme_remote_version(database, scheme_id, session=session)
    except Exception:
        return False
    loc = meta.get('remote_snapshot') or {}
    keys = ('records', 'locus_count', 'last_updated', 'last_added')
    compared = False
    for k in keys:
        r = remote.get(k)
        l = loc.get(k)
        if r is None:
            continue
        if l is None:
            return False
        compared = True
        if l != r:
            return False
    return compared


def fetch_public_scheme_catalog(session=None):
    root_resources = fetch_json(f'{REST_BASE}/db', session=session)
    databases = []
    for group in root_resources:
        databases.extend(group.get('databases', []))

    seqdef_databases = []
    for db in databases:
        href = db.get('href')
        name = db.get('name', '')
        if href and name.endswith('_seqdef'):
            seqdef_databases.append(db)

    catalog = []
    progress = tqdm(
        seqdef_databases,
        desc='Discovering schemes from API',
        unit='database',
        leave=True
    )
    for db in progress:
        db_name = db.get('name')
        progress.set_postfix_str(db_name)
        schemes_url = f'{REST_BASE}/db/{db_name}/schemes?with_pk=1'
        schemes = fetch_json(schemes_url, session=session).get('schemes', [])
        for entry in schemes:
            scheme_url = entry.get('scheme')
            if not scheme_url:
                continue
            try:
                detail = fetch_json(scheme_url, session=session)
            except Exception:
                continue
            loci_urls = detail.get('loci') or []
            if not loci_urls:
                loci_urls = fetch_json(f'{scheme_url}/loci', session=session).get('loci', [])
            description = (detail.get('description') or entry.get('description') or '')
            if not loci_urls:
                continue

            scheme_id_int = int(scheme_url.rstrip('/').split('/')[-1])
            codename = ensure_scheme_codename(db_name, scheme_id_int, description)
            species = db_name
            if species.startswith('pubmlst_'):
                species = species[len('pubmlst_'):]
            if species.endswith('_seqdef'):
                species = species[:-len('_seqdef')]

            catalog.append({
                'codename': codename,
                'description': description.strip(),
                'database': db_name,
                'species': species,
                'scheme_id': scheme_id_int,
                'profiles_csv': f'{scheme_url}/profiles_csv',
                'loci': loci_urls,
                'remote_version': remote_version_from_scheme_detail(detail),
            })
        progress.set_postfix_str(f'{db_name} | total schemes: {len(catalog)}')
    return catalog


def parse_direct_scheme_selector(selector):
    """
    Resolve a scheme selector to (database, scheme_id):

    * ``database:scheme_id`` — e.g. ``pubmlst_cdifficile_seqdef:1``
    * Stable codename suffix — ``database_schemeid`` with numeric id at the end,
      e.g. ``pubmlst_cdifficile_seqdef_1`` (same layout as the on-disk folder).
    """
    s = (selector or '').strip()
    if not s:
        return None, None
    if ':' in s:
        db, sid = s.rsplit(':', 1)
        db = db.strip()
        sid = sid.strip()
        if not db or not sid.isdigit():
            return None, None
        return db, int(sid)
    m = re.fullmatch(r'(.+)_(\d+)', s)
    if not m:
        return None, None
    db = m.group(1).strip()
    if not db:
        return None, None
    return db, int(m.group(2))


def fetch_scheme_item_direct(database, scheme_id, session=None):
    """Load one scheme from the API without scanning all databases."""
    scheme_url = f'{REST_BASE}/db/{database}/schemes/{scheme_id}'
    try:
        detail = fetch_json(scheme_url, session=session)
    except Exception as exc:
        raise RuntimeError(
            f"Scheme not found or not accessible: {database}:{scheme_id}. "
            'Check the database name and scheme id.'
        ) from exc

    loci_urls = detail.get('loci') or []
    if not loci_urls:
        try:
            loci_payload = fetch_json(f'{scheme_url}/loci', session=session)
            loci_urls = loci_payload.get('loci', [])
        except Exception as exc:
            raise RuntimeError(
                f"Could not load loci for scheme {database}:{scheme_id}."
            ) from exc
    if not loci_urls:
        raise RuntimeError(
            f"Scheme {database}:{scheme_id} has no loci in the API response."
        )

    description = (detail.get('description') or '').strip()
    codename = ensure_scheme_codename(database, scheme_id, description)
    species = database
    if species.startswith('pubmlst_'):
        species = species[len('pubmlst_'):]
    if species.endswith('_seqdef'):
        species = species[:-len('_seqdef')]

    return {
        'codename': codename,
        'description': description,
        'database': database,
        'species': species,
        'scheme_id': scheme_id,
        'profiles_csv': f'{scheme_url}/profiles_csv',
        'loci': loci_urls,
        'remote_version': remote_version_from_scheme_detail(detail),
    }


def build_scheme_catalog_direct(selectors, session=None):
    """
    Build a minimal catalog from selectors (``database:scheme_id`` or
    ``database_<id>`` codename). No GET /db crawl and no per-database listing.
    """
    catalog = []
    seen_keys = set()
    for raw in selectors:
        sel = raw.strip()
        if not sel:
            continue
        db, sid = parse_direct_scheme_selector(sel)
        if db is None:
            raise RuntimeError(
                f'Invalid selector {sel!r}. Use database:scheme_id or the stable codename '
                '(e.g. pubmlst_cdifficile_seqdef:1 or pubmlst_cdifficile_seqdef_1).'
            )
        key = (db, sid)
        if key in seen_keys:
            continue
        seen_keys.add(key)
        item = fetch_scheme_item_direct(db, sid, session=session)
        catalog.append(item)
    return catalog


def save_scheme_catalog(catalog):
    pathdb.mkdir(parents=True, exist_ok=True)
    (pathdb / SCHEME_CATALOG_FILE).write_text(
        json.dumps(catalog, indent=2),
        encoding='utf-8'
    )
    meta = {
        'updated_at': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        'count': len(catalog),
    }
    (pathdb / SCHEME_CATALOG_META_FILE).write_text(
        json.dumps(meta, indent=2),
        encoding='utf-8'
    )


def load_scheme_catalog():
    catalog_path = pathdb / SCHEME_CATALOG_FILE
    if catalog_path.is_file():
        return json.loads(catalog_path.read_text(encoding='utf-8'))
    bundled = _package_bundle_dir() / SCHEME_CATALOG_FILE
    if bundled.is_file():
        return json.loads(bundled.read_text(encoding='utf-8'))
    return None


def load_scheme_catalog_meta():
    meta_path = pathdb / SCHEME_CATALOG_META_FILE
    if meta_path.is_file():
        try:
            return json.loads(meta_path.read_text(encoding='utf-8'))
        except Exception:
            return None
    bundled = _package_bundle_dir() / SCHEME_CATALOG_META_FILE
    if bundled.is_file():
        try:
            return json.loads(bundled.read_text(encoding='utf-8'))
        except Exception:
            return None
    return None


def scheme_catalog_is_user_cache():
    """True if the active catalog JSON is the copy under pathdb (not the package bundle)."""
    return (pathdb / SCHEME_CATALOG_FILE).is_file()


def get_remote_scheme_catalog(session=None, force_refresh=False):
    """
    Return (catalog, from_cache, used_live_api).

    * ``used_live_api`` is True when ``fetch_public_scheme_catalog`` ran in this call.
    * If ``force_refresh`` is True, the live API is always used (never cache/bundle only).
    """
    pathdb.mkdir(parents=True, exist_ok=True)
    if force_refresh:
        catalog = fetch_public_scheme_catalog(session=session)
        save_scheme_catalog(catalog)
        return catalog, False, True
    cached = load_scheme_catalog()
    if cached:
        return cached, True, False
    catalog = fetch_public_scheme_catalog(session=session)
    save_scheme_catalog(catalog)
    return catalog, False, True


def print_remote_scheme_catalog(catalog, from_cache=False, meta=None):
    suffix = ''
    if from_cache:
        if scheme_catalog_is_user_cache():
            if meta and meta.get('updated_at'):
                suffix = f" (cached, updated {meta['updated_at']})"
            else:
                suffix = ' (cached)'
        else:
            if meta and meta.get('updated_at'):
                suffix = f" (bundled with package, snapshot {meta['updated_at']})"
            else:
                suffix = ' (bundled with package)'
    else:
        if meta and meta.get('updated_at'):
            suffix = f" (refreshed {meta['updated_at']})"
    print(f'Total remote schemes: {len(catalog)}{suffix}\n')
    for i, item in enumerate(catalog, start=1):
        print(
            f"({i}) {item['codename']} | {item['database']}:{item['scheme_id']} | "
            f"{item.get('species', '')} | {item['description']}"
        )


def download_scheme_data(items):
    """
    Download one scheme. Never raises: returns
    (ok, codename, loci_names, description, failed_item, error_message).
    On failure, removes the scheme directory to avoid a half-written tree.
    """
    item, session = items
    codename = (item.get('codename') or '').strip()
    outdir = pathdb / 'schemes' / codename
    db = item.get('database')
    sid = item.get('scheme_id')
    label = f'{db}:{sid}' if db is not None and sid is not None else repr(item.get('codename'))

    if not codename:
        msg = f'Internal error: empty scheme codename for {label}'
        logger.error(msg)
        return (False, None, None, None, item, msg)

    progress = None
    try:
        outdir.mkdir(exist_ok=True, parents=True)
        expected_files = 1 + len(item.get('loci') or [])
        bytes_downloaded = 0
        progress = tqdm(
            total=expected_files,
            desc=f'Downloading files for {codename}',
            unit='file',
            leave=False,
        )

        profiles = fetch_text(
            item['profiles_csv'],
            content_type='text/tab-separated-values',
            session=session,
        )
        (outdir / f'{codename}.txt').write_text(profiles, encoding='utf-8')
        bytes_downloaded += _text_nbytes(profiles)
        progress.update(1)
        progress.set_postfix_str(f'{bytes_downloaded / (1024 * 1024):.2f} MiB')

        loci_names = []
        loci_without_fasta = 0
        for locus_url in item['loci']:
            locus_record = fetch_json(locus_url, session=session)
            locus_name = locus_record.get('id') or infer_locus_name_from_url(locus_url)
            alleles_fasta = locus_record.get('alleles_fasta')
            if not alleles_fasta:
                loci_without_fasta += 1
                logger.warning(
                    'No alleles_fasta link for locus %r in scheme %s (%s)',
                    locus_name,
                    label,
                    codename,
                )
                progress.update(1)
                progress.set_postfix_str(f'{bytes_downloaded / (1024 * 1024):.2f} MiB')
                continue
            fasta_data = fetch_text(alleles_fasta, content_type='text/x-fasta', session=session)
            (outdir / f'{locus_name}.tfa').write_text(fasta_data, encoding='utf-8')
            bytes_downloaded += _text_nbytes(fasta_data)
            progress.update(1)
            progress.set_postfix_str(f'{bytes_downloaded / (1024 * 1024):.2f} MiB')
            loci_names.append(locus_name)
        if loci_without_fasta:
            msg = (
                f'Incomplete allele FASTA download for scheme {label} ({codename}). '
                f'{loci_without_fasta} of {len(item.get("loci") or [])} locus/loci had no '
                'alleles_fasta URL (often fixed with PubMLST OAuth).'
            )
            raise RuntimeError(msg)
        if not loci_names:
            msg = (
                f'No allele FASTA could be downloaded for scheme {label} ({codename}). '
                f'{loci_without_fasta} locus/loci had no alleles_fasta URL (often fixed with PubMLST OAuth).'
            )
            raise RuntimeError(msg)

        snap = item.get('remote_version') or fetch_scheme_remote_version(
            item['database'], item['scheme_id'], session=session
        )
        save_local_scheme_meta(
            codename, item['database'], item['scheme_id'], item['description'], snap
        )

        return (True, codename, loci_names, item['description'], None, None)
    except Exception as exc:
        msg = str(exc)
        logger.error('Scheme download failed for %s (%s): %s', label, codename, msg)
        if outdir.is_dir():
            shutil.rmtree(outdir, ignore_errors=True)
        return (False, None, None, None, item, msg)
    finally:
        if progress is not None:
            progress.close()


def parse_update_mlst_selectors(raw_selectors):
    """
    Validate selectors passed to ``--update-mlst``.

    Returns ``None`` when the user requests the full catalog (token ``ALL``),
    otherwise a non-empty list of selector strings.

    Raises:
        ValueError: with a message suitable for ``argparse.ArgumentParser.error``.
    """
    if not raw_selectors:
        raise ValueError(
            '--update-mlst requires a value. Use ALL for the full catalog (very slow) '
            'or comma-separated database:scheme_id / stable codenames (see --scheme-list).'
        )
    if len(raw_selectors) == 1 and raw_selectors[0].upper() == 'ALL':
        return None
    if any(s.upper() == 'ALL' for s in raw_selectors):
        raise ValueError('ALL cannot be combined with other --update-mlst values.')
    return raw_selectors


def update_mlstdb(threads):
    return update_mlstdb_selected(threads=threads, selectors=None)


def _load_or_build_scheme_catalog(session=None):
    """
    For full DB updates: reuse scheme_catalog.json from pathdb or package bundle
    if present; otherwise fetch from the API and save under pathdb.
    Returns (catalog, used_cache).
    """
    cached = load_scheme_catalog()
    if cached:
        return cached, True
    catalog = fetch_public_scheme_catalog(session=session)
    save_scheme_catalog(catalog)
    return catalog, False


def update_mlstdb_selected(threads, selectors=None, force_redownload_schemes=False):
    pathdb.mkdir(exist_ok=True, parents=True)
    session = _build_authenticated_session()
    if session is not None:
        logger.info('Using authenticated PubMLST API access (OAuth1)')
        if threads > 1:
            logger.info('OAuth session is used in single-thread mode to avoid concurrent token/session issues.')
            threads = 1
    else:
        logger.warning('Using anonymous PubMLST API access; post-2024 records may be unavailable.')
    if selectors:
        logger.info(
            'Resolving selector(s) directly as database:scheme_id (no full catalog scan)'
        )
        scheme_catalog = build_scheme_catalog_direct(selectors, session=session)
        logger.info('Resolved %s scheme(s) to download', len(scheme_catalog))
        logger.info('Selectors: %s', ', '.join(selectors))
    else:
        scheme_catalog, catalog_from_cache = _load_or_build_scheme_catalog(session=session)
        if catalog_from_cache:
            src = 'user cache' if scheme_catalog_is_user_cache() else 'package bundle'
            logger.info(
                'Using PubMLST scheme catalog from %s (%s schemes). '
                'To rebuild from the API: fastmlst --scheme-list-update '
                '(for downloads use --update-mlst … or ALL).',
                src,
                len(scheme_catalog),
            )
        else:
            logger.info(
                'Built scheme catalog from BIGSdb REST API (%s schemes)',
                len(scheme_catalog),
            )
    scheme_number = defaultdict(list)
    scheme_list = {}
    if selectors and (pathdb / 'scheme_number.pkl').is_file():
        try:
            scheme_number.update(load_obj(str(pathdb / 'scheme_number.pkl')))
        except Exception:
            pass
    if selectors and (pathdb / SCHEME_LIST_FILE).is_file():
        try:
            scheme_list.update(load_obj(str(pathdb / SCHEME_LIST_FILE)))
        except Exception:
            pass

    to_download = []
    if force_redownload_schemes:
        to_download = list(scheme_catalog)
    else:
        for item in scheme_catalog:
            cn = item['codename']
            if local_scheme_matches_remote(
                cn, item['database'], item['scheme_id'], session=session
            ):
                logger.info(
                    'Skipping %s (%s:%s) — local copy matches remote metadata',
                    cn, item['database'], item['scheme_id'],
                )
                scheme_number[cn] = local_scheme_locus_stems(cn)
                lm = load_local_scheme_meta(cn)
                scheme_list[cn] = (
                    (lm.get('description') if lm else None)
                    or item.get('description', '')
                )
            else:
                to_download.append(item)

    download_failures = []
    if to_download:
        logger.info(
            'Downloading %s scheme(s) (%s unchanged)',
            len(to_download),
            len(scheme_catalog) - len(to_download),
        )
        with ThreadPool(threads) as t:
            download_items = [(item, session) for item in to_download]
            for ok, codename, loci_names, description, failed_item, err_msg in tqdm(
                t.imap(download_scheme_data, download_items),
                total=len(to_download),
                desc='Downloading Schemes using {} threads'.format(threads),
                unit='Schemes',
                leave=True,
            ):
                if ok:
                    scheme_number[codename] = loci_names
                    scheme_list[codename] = description
                else:
                    download_failures.append((failed_item, err_msg))
                    failed_cn = (failed_item.get('codename') or '').strip()
                    if failed_cn:
                        scheme_number.pop(failed_cn, None)
                        scheme_list.pop(failed_cn, None)
        if download_failures:
            logger.warning(
                '%d scheme download(s) failed; continuing with %d scheme(s) on disk',
                len(download_failures),
                len(scheme_number),
            )
    else:
        logger.info(
            'All %s scheme(s) are up to date with the remote API; no downloads',
            len(scheme_catalog),
        )

    if not scheme_number:
        if download_failures:
            detail = '; '.join(
                f"{it.get('database')}:{it.get('scheme_id')}: {err}"
                for it, err in download_failures
            )
            raise RuntimeError(
                'No MLST schemes were available after update. All downloads failed. ' + detail
            )
        raise RuntimeError('No MLST schemes were available after update (nothing on disk).')

    # Backwards-compatibility marker file used by safety checks.
    (pathdb / DBASES_COMPAT_FILE).write_text(
        '<resources source="https://rest.pubmlst.org/db"/>',
        encoding='utf-8'
    )
    logger.info('Schemes were downloaded')
    schemes_root = (pathdb / 'schemes').resolve()
    allfasta = []
    save_obj(dict(scheme_number), str(pathdb) + '/scheme_number.pkl')
    save_obj(dict(scheme_list), str(pathdb) + '/' + SCHEME_LIST_FILE)
    logger.info('Schemes object was created in {}'.format(str(pathdb) +
                                                          '/scheme_number.pkl')
                )
    skipped_mlst_fasta = set()
    for fasta in sorted(schemes_root.rglob('*.tfa')):
        try:
            rel = fasta.resolve().relative_to(schemes_root)
        except ValueError:
            continue
        if len(rel.parts) >= 2:
            scheme = rel.parts[0]
        else:
            scheme = rel.stem
        if scheme_excluded_from_mlst_concat(scheme):
            skipped_mlst_fasta.add(scheme)
            continue
        for record in SeqIO.parse(str(fasta), 'fasta'):
            record.id = '{}.{}'.format(scheme, record.id)
            record.description = ''
            allfasta.append(record)
    if skipped_mlst_fasta:
        logger.info(
            'Omitted from mlst.fasta merge: %s',
            ', '.join(sorted(skipped_mlst_fasta)),
        )
    if not allfasta:
        raise RuntimeError(
            'No allele sequences were found under schemes/**/*.tfa; cannot build mlst.fasta. '
            'If you previously had an empty scheme codename, remove stray *.tfa under schemes/ '
            'and run the update again. '
            'Note: pubmlst_rmlst_seqdef schemes are never merged into mlst.fasta; install at least '
            'one standard MLST scheme alongside rMLST if you need a BLAST allele database.'
        )
    outfna = 'mlst.fasta'
    SeqIO.write(allfasta, str(pathdb) + '/' + outfna, 'fasta')
    blastdb_cmd = [
        'makeblastdb',
        '-hash_index',
        '-in', str(pathdb / outfna),
        '-dbtype', 'nucl',
        '-title', f'PubMLST_{date.today().strftime("%d%m%y")}',
        '-parse_seqids'
    ]
    try:
        subprocess.run(
            blastdb_cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
    except subprocess.CalledProcessError as exc:
        logger.error('makeblastdb failed while building PubMLST database')
        if exc.stderr:
            logger.error(exc.stderr.strip())
        raise
    logger.info('BLASTdb was created using pubmlst data')
    logger.info('Update PubMLST Complete')

def _count_fasta_records(fasta_path):
    """Count sequences in a FASTA file (lines starting with '>')."""
    with fasta_path.open('rb') as fh:
        return sum(1 for line in fh if line.startswith(b'>'))


def scheme_type_from_description(description):
    """
    Short label for PubMLST scheme descriptions (MLST, cgMLST, etc.).
    """
    d = (description or '').strip().lower()
    if not d:
        return 'unknown'
    if 'cgmlst' in d or 'cg-mlst' in d:
        return 'cgMLST'
    if 'wgmlst' in d or 'whole genome' in d and 'mlst' in d:
        return 'wgMLST'
    if 'mlst' in d:
        return 'MLST'
    if 'typing' in d:
        return 'typing'
    return 'other'


def _count_st_rows_in_profile(profile_path):
    """
    Number of ST (profile) rows in a PubMLST profiles TSV: non-empty data lines
    after the header.
    """
    if not profile_path.is_file():
        return None
    text = profile_path.read_text(encoding='utf-8', errors='replace')
    lines = [
        ln for ln in text.splitlines()
        if ln.strip() and not ln.lstrip().startswith('#')
    ]
    if not lines:
        return 0
    return max(0, len(lines) - 1)


def _stats_description_one_line(description):
    """Collapse whitespace and strip; safe for a single pipe-separated output line."""
    text = ' '.join((description or '').split())
    return text.replace('|', '/')


def print_installed_scheme_stats():
    """
    Print one line per scheme under pathdb/schemes: API description, inferred type,
    ST count, total alleles, and per-locus allele counts.
    """
    schemes_root = pathdb / 'schemes'
    if not schemes_root.is_dir():
        print(f'No schemes directory: {schemes_root}')
        return
    subdirs = sorted(p for p in schemes_root.iterdir() if p.is_dir())
    if not subdirs:
        print(f'No scheme folders under {schemes_root}')
        return
    scheme_descriptions = {}
    pkl_path = pathdb / SCHEME_LIST_FILE
    if pkl_path.is_file():
        try:
            scheme_descriptions.update(load_obj(str(pkl_path)))
        except Exception:
            pass
    for scheme_dir in subdirs:
        cn = scheme_dir.name
        profile = scheme_dir / f'{cn}.txt'
        st_n = _count_st_rows_in_profile(profile)
        tfa_files = sorted(scheme_dir.glob('*.tfa'))
        per_locus = []
        total_alleles = 0
        for tfa in tfa_files:
            n = _count_fasta_records(tfa)
            total_alleles += n
            per_locus.append(f'{tfa.stem}={n}')
        st_str = str(st_n) if st_n is not None else '?'
        locus_part = ', '.join(per_locus) if per_locus else '(no .tfa)'
        desc = (scheme_descriptions.get(cn) or '').strip()
        if not desc:
            meta = load_local_scheme_meta(cn)
            if meta:
                desc = (meta.get('description') or '').strip()
        stype = scheme_type_from_description(desc)
        desc_line = _stats_description_one_line(desc)
        if not desc_line:
            desc_line = '(none)'
        print(
            f'{cn} | type={stype} | description={desc_line} | STs={st_str} | '
            f'alleles={total_alleles} | per_locus: {locus_part}'
        )


def show_scheme_list():
    catalog = load_scheme_catalog()
    if catalog is not None:
        print_remote_scheme_catalog(
            catalog,
            from_cache=True,
            meta=load_scheme_catalog_meta(),
        )
        return
    scheme_list = pathdb / SCHEME_LIST_FILE
    if not scheme_list.is_file():
        from sys import exit
        logger.error('There is no cached scheme list, please update the database')
        exit()
    species = load_obj(str(scheme_list))
    print(f'There are {len(species)} schemes (A round of applause to @keithajolley! (Jolley, et al., 2018)):\n')
    i = 1
    for sch, specie in species.items():
        print(f'({i}) {sch}: {specie.strip()}')
        i += 1
