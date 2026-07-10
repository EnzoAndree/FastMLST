import json
import logging
import os
import re
import shutil
import stat
import tempfile
import time
import uuid
from collections import defaultdict
from contextlib import contextmanager
from datetime import date, datetime, timezone
from io import StringIO
from multiprocessing.pool import ThreadPool
from pathlib import Path
from pickle import dump
from pickle import load
from urllib.parse import unquote, urlparse
from urllib.request import Request, urlopen

import subprocess
from Bio import SeqIO
from tqdm import tqdm  # pip3 install tqdm

from fastmlst.pubmlst_client import PubMLSTAuthenticationError
from fastmlst.pubmlst_client import PubMLSTClient
from fastmlst.pubmlst_client import PubMLSTHTTPError

logger = logging.getLogger('update_mlst')
REST_BASE = 'https://rest.pubmlst.org'
DBASES_COMPAT_FILE = 'dbases.xml'
SCHEME_LIST_FILE = 'scheme_list.pkl'
SCHEME_CATALOG_FILE = 'scheme_catalog.json'
SCHEME_CATALOG_META_FILE = 'scheme_catalog_meta.json'
SCHEME_LOCAL_META_FILE = '.fastmlst_scheme_meta.json'
UPDATE_LOCK_SUFFIX = '.fastmlst-update.lock'
LARGE_SCHEME_LOCUS_THRESHOLD = 100
BLAST_INDEX_SUFFIXES = (
    '.nal', '.ndb', '.nhr', '.nin', '.njs', '.nog', '.nos', '.not', '.nsq',
    '.ntf', '.nto',
)
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


def _validate_pubmlst_database_name(database):
    database = str(database or '').strip()
    if not re.fullmatch(r'[A-Za-z0-9_]+', database):
        raise ValueError(f'Invalid PubMLST database name: {database!r}')
    return database

# NEW: Define a function to override the default database location
def set_pathdb(custom_path):
    """
    Override the default PubMLST database directory.

    Parameters:
        custom_path (str): Custom path for the MLST database directory.
    """
    global pathdb
    pathdb = Path(custom_path).expanduser()


@contextmanager
def _using_database_path(database_path):
    """Temporarily redirect module-level database helpers to a staging tree."""
    global pathdb
    previous = pathdb
    pathdb = Path(database_path)
    try:
        yield pathdb
    finally:
        pathdb = previous


def _validate_database_destination(database_path):
    resolved = Path(database_path).expanduser().resolve()
    protected_paths = (Path.home().resolve(), Path.cwd().resolve())
    if resolved == Path('/').resolve():
        raise RuntimeError(f'Refusing to use unsafe database path: {resolved}')
    for protected in protected_paths:
        try:
            protected.relative_to(resolved)
        except ValueError:
            continue
        raise RuntimeError(
            f'Refusing to use database path containing a protected directory: {resolved}'
        )
    return resolved


def _make_lock_file_readable(lock_path):
    """Ensure every database reader can open the persistent coordination lock."""
    try:
        mode = stat.S_IMODE(lock_path.stat().st_mode)
        required_mode = mode | stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH
        if required_mode != mode:
            lock_path.chmod(required_mode)
    except OSError as exc:
        raise RuntimeError(
            f'Could not make the FastMLST lock file readable: {lock_path}'
        ) from exc


@contextmanager
def _database_update_lock(database_path):
    """Prevent concurrent writers from publishing competing database trees."""
    database_path = _validate_database_destination(database_path)
    database_path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = database_path.parent / f'.{database_path.name}{UPDATE_LOCK_SUFFIX}'
    handle = lock_path.open('a+')
    lock_backend = None
    try:
        _make_lock_file_readable(lock_path)
        try:
            import fcntl
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            lock_backend = ('fcntl', fcntl)
        except ImportError:
            try:
                import msvcrt
                handle.seek(0, os.SEEK_END)
                if handle.tell() == 0:
                    handle.write('\0')
                    handle.flush()
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
                lock_backend = ('msvcrt', msvcrt)
            except (ImportError, OSError) as exc:
                raise RuntimeError(
                    f'Could not acquire the FastMLST update lock for {database_path}'
                ) from exc
        except BlockingIOError as exc:
            raise RuntimeError(
                f'Another FastMLST update is already using {database_path}'
            ) from exc
        handle.seek(0)
        handle.truncate()
        handle.write(f'pid={os.getpid()} started={datetime.now(timezone.utc).isoformat()}\n')
        handle.flush()
        yield
    finally:
        try:
            if lock_backend and lock_backend[0] == 'fcntl':
                lock_backend[1].flock(handle.fileno(), lock_backend[1].LOCK_UN)
            elif lock_backend and lock_backend[0] == 'msvcrt':
                handle.seek(0)
                lock_backend[1].locking(handle.fileno(), lock_backend[1].LK_UNLCK, 1)
        except OSError:
            pass
        handle.close()


def acquire_database_read_lock(database_root=None):
    """Acquire a shared reader lock; caller must keep and eventually close the handle."""
    database_root = Path(database_root or pathdb).expanduser().resolve()
    lock_path = database_root.parent / f'.{database_root.name}{UPDATE_LOCK_SUFFIX}'

    # Readers commonly use a centrally installed, read-only database. Open an
    # existing lock file read-only so those installations remain usable.
    try:
        handle = lock_path.open('r')
    except FileNotFoundError:
        try:
            database_root.parent.mkdir(parents=True, exist_ok=True)
            handle = lock_path.open('a+')
        except OSError as exc:
            raise RuntimeError(
                f'Could not initialize the FastMLST read lock for {database_root}. '
                f'For a shared read-only database, ask its administrator to create '
                f'a readable sibling lock file at {lock_path}.'
            ) from exc
        try:
            _make_lock_file_readable(lock_path)
        except Exception:
            handle.close()
            raise
    except OSError as exc:
        raise RuntimeError(
            f'Could not open the FastMLST read lock for {database_root}'
        ) from exc

    try:
        try:
            import fcntl
            fcntl.flock(handle.fileno(), fcntl.LOCK_SH | fcntl.LOCK_NB)
            handle._fastmlst_lock_backend = ('fcntl', fcntl)
        except ImportError:
            import msvcrt
            handle.seek(0, os.SEEK_END)
            if handle.tell() == 0:
                if not handle.writable():
                    raise OSError('The read-only lock file is empty')
                handle.write('\0')
                handle.flush()
            handle.seek(0)
            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            handle._fastmlst_lock_backend = ('msvcrt', msvcrt)
    except (BlockingIOError, OSError) as exc:
        handle.close()
        raise RuntimeError(
            f'The FastMLST database is being updated: {database_root}'
        ) from exc
    return handle


def release_database_read_lock(handle):
    if handle is None or handle.closed:
        return
    backend = getattr(handle, '_fastmlst_lock_backend', None)
    try:
        if backend and backend[0] == 'fcntl':
            backend[1].flock(handle.fileno(), backend[1].LOCK_UN)
        elif backend and backend[0] == 'msvcrt':
            handle.seek(0)
            backend[1].locking(handle.fileno(), backend[1].LK_UNLCK, 1)
    finally:
        handle.close()


@contextmanager
def database_read_lock(database_root=None):
    handle = acquire_database_read_lock(database_root)
    try:
        yield
    finally:
        release_database_read_lock(handle)


def _link_or_copy(source, destination):
    """Use cheap copy-on-write hard links, falling back to a real copy."""
    try:
        os.link(str(source), str(destination))
        return destination
    except OSError:
        return shutil.copy2(str(source), str(destination))


def _prepare_database_staging(database_path, clone_existing):
    database_path = _validate_database_destination(database_path)
    stage = Path(tempfile.mkdtemp(
        prefix=f'.{database_path.name}.staging-', dir=str(database_path.parent)
    ))
    if clone_existing and database_path.is_dir():
        stage.rmdir()
        shutil.copytree(
            str(database_path),
            str(stage),
            symlinks=True,
            copy_function=_link_or_copy,
        )
    elif database_path.is_dir():
        for filename in (SCHEME_CATALOG_FILE, SCHEME_CATALOG_META_FILE):
            source = database_path / filename
            if source.is_file():
                shutil.copy2(str(source), str(stage / filename))
    return stage


def _recover_interrupted_database_swap(database_path):
    """Restore a backup left between the two directory renames of a swap."""
    database_path = _validate_database_destination(database_path)
    backups = sorted(
        database_path.parent.glob(f'.{database_path.name}.backup-*'),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    if backups:
        backup = backups[0]
        target_usable = blast_database_is_ready(database_path, deep=True)
        backup_has_data = backup.is_dir() and any(backup.iterdir())
        if not database_path.exists() or (not target_usable and backup_has_data):
            if database_path.exists():
                shutil.rmtree(str(database_path), ignore_errors=False)
            os.replace(str(backup), str(database_path))
            backups = backups[1:]
    for backup in backups:
        shutil.rmtree(str(backup), ignore_errors=True)
    for stale_stage in database_path.parent.glob(f'.{database_path.name}.staging-*'):
        shutil.rmtree(str(stale_stage), ignore_errors=True)


def _publish_database_staging(database_path, stage):
    """Publish a verified tree and roll back if the second rename fails."""
    database_path = _validate_database_destination(database_path)
    stage = Path(stage).resolve()
    if stage.parent != database_path.parent or not stage.is_dir():
        raise RuntimeError('Database staging directory is invalid or on another filesystem.')
    backup = database_path.parent / (
        f'.{database_path.name}.backup-{uuid.uuid4().hex}'
    )
    moved_old = False
    try:
        if database_path.exists():
            os.replace(str(database_path), str(backup))
            moved_old = True
        os.replace(str(stage), str(database_path))
    except Exception:
        if moved_old and backup.exists() and not database_path.exists():
            os.replace(str(backup), str(database_path))
        raise
    if backup.exists():
        shutil.rmtree(str(backup), ignore_errors=True)


def recover_database_if_needed(database_root=None):
    """Recover the active database after a process died between directory renames."""
    database_root = Path(database_root or pathdb).expanduser().resolve()
    backups = list(database_root.parent.glob(f'.{database_root.name}.backup-*'))
    if not backups:
        return False
    database_root = _validate_database_destination(database_root)
    with _database_update_lock(database_root):
        _recover_interrupted_database_swap(database_root)
    return True

# Define the path to the .cache directory in the user's home directory
home_dir = Path.home()
cache_dir = home_dir / '.cache' / 'fastmlst'

# Ensure the cache directory exists
cache_dir.mkdir(parents=True, exist_ok=True)

# Default pathdb points to the pubmlst folder in the default cache
pathdb = cache_dir / 'pubmlst'

necessary_file = ['mlst.fasta.nhr', 'mlst.fasta.nsq', 'mlst.fasta.nin']
_oauth_credentials = None


def _atomic_write_text(destination, text, encoding='utf-8'):
    """Write a text file without exposing a partially-written destination."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f'.{destination.name}.', suffix='.tmp', dir=str(destination.parent)
    )
    try:
        with os.fdopen(fd, 'w', encoding=encoding) as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_name, str(destination))
    except Exception:
        try:
            os.unlink(tmp_name)
        except FileNotFoundError:
            pass
        raise


def _atomic_save_obj(obj, destination):
    """Pickle an object to a temporary file and atomically publish it."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f'.{destination.name}.', suffix='.tmp', dir=str(destination.parent)
    )
    try:
        with os.fdopen(fd, 'wb') as handle:
            dump(obj, handle, 2)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_name, str(destination))
    except Exception:
        try:
            os.unlink(tmp_name)
        except FileNotFoundError:
            pass
        raise


def save_obj(obj, name):
    _atomic_save_obj(obj, name)


def load_obj(name):
    with open(name, 'rb') as f:
        return load(f)


def _ensure_config_dir():
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    os.chmod(str(CONFIG_DIR), 0o700)


def save_pubmlst_client_credentials(client_id, client_secret, auth_db=None):
    if not client_id or not client_secret:
        raise ValueError('client_id and client_secret are required.')
    _ensure_config_dir()
    payload = {
        'client_id': client_id,
        'client_secret': client_secret,
        'auth_db': _validate_pubmlst_database_name(auth_db or AUTH_TEST_DB),
    }
    _atomic_write_text(CONFIG_FILE, json.dumps(payload, indent=2), encoding='utf-8')
    os.chmod(str(CONFIG_FILE), 0o600)


def load_pubmlst_client_credentials():
    if not CONFIG_FILE.is_file():
        return None, None
    try:
        payload = json.loads(CONFIG_FILE.read_text(encoding='utf-8'))
    except Exception:
        return None, None
    return payload.get('client_id'), payload.get('client_secret')


def load_pubmlst_auth_database():
    if not CONFIG_FILE.is_file():
        return AUTH_TEST_DB
    try:
        payload = json.loads(CONFIG_FILE.read_text(encoding='utf-8'))
    except Exception:
        return AUTH_TEST_DB
    return payload.get('auth_db') or AUTH_TEST_DB


def set_oauth_credentials(client_id=None, client_secret=None,
                          access_token=None, access_secret=None,
                          verifier=None, auth_db=None):
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
        _oauth_credentials = {
            'client_id': client_id,
            'client_secret': client_secret,
            'auth_db': _validate_pubmlst_database_name(
                auth_db or load_pubmlst_auth_database()
            ),
        }
        if access_token and access_secret:
            _oauth_credentials['access_token'] = access_token
            _oauth_credentials['access_secret'] = access_secret
        if verifier:
            _oauth_credentials['verifier'] = verifier
    else:
        _oauth_credentials = None


def configure_pubmlst_auth_from_environment(require_auth=False, **overrides):
    """
    Load OAuth inputs from explicit overrides, environment variables, or the
    protected configuration file. Credentials supplied through the environment
    are deliberately not persisted.
    """
    client_id = overrides.get('client_id') or os.environ.get(
        'FASTMLST_PUBMLST_CLIENT_ID'
    )
    client_secret = overrides.get('client_secret') or os.environ.get(
        'FASTMLST_PUBMLST_CLIENT_SECRET'
    )
    access_token = overrides.get('access_token') or os.environ.get(
        'FASTMLST_PUBMLST_ACCESS_TOKEN'
    )
    access_secret = overrides.get('access_secret') or os.environ.get(
        'FASTMLST_PUBMLST_ACCESS_SECRET'
    )
    verifier = overrides.get('verifier') or os.environ.get(
        'FASTMLST_PUBMLST_VERIFIER'
    )
    auth_db = overrides.get('auth_db') or os.environ.get(
        'FASTMLST_PUBMLST_AUTH_DB'
    )
    set_oauth_credentials(
        client_id=client_id,
        client_secret=client_secret,
        access_token=access_token,
        access_secret=access_secret,
        verifier=verifier,
        auth_db=auth_db,
    )
    if require_auth and not _oauth_credentials:
        raise RuntimeError(
            'Authenticated PubMLST access is required. Run fastmlst '
            '--pubmlst-connect with your client ID and secret first.'
        )
    return _oauth_credentials is not None


def _oauth_cache_path():
    _ensure_config_dir()
    return TOKENS_FILE


def _load_cached_access_token(client_id, auth_db):
    cache = _oauth_cache_path()
    if not cache.is_file():
        return None, None
    try:
        payload = json.loads(cache.read_text(encoding='utf-8'))
    except Exception:
        return None, None
    record = payload.get(client_id, {})
    # Access tokens are granted for a specific BIGSdb resource. Legacy cache
    # records did not store that resource, so they must be authorized again.
    if record.get('auth_db') != auth_db:
        return None, None
    return record.get('access_token'), record.get('access_secret')


def _save_cached_access_token(client_id, access_token, access_secret, auth_db):
    cache = _oauth_cache_path()
    payload = {}
    if cache.is_file():
        try:
            payload = json.loads(cache.read_text(encoding='utf-8'))
        except Exception:
            payload = {}
    payload[client_id] = {
        'access_token': access_token,
        'access_secret': access_secret,
        'auth_db': _validate_pubmlst_database_name(auth_db),
    }
    _atomic_write_text(cache, json.dumps(payload, indent=2), encoding='utf-8')
    os.chmod(str(cache), 0o600)


def _delete_cached_access_token(client_id):
    cache = _oauth_cache_path()
    if not cache.is_file():
        return
    try:
        payload = json.loads(cache.read_text(encoding='utf-8'))
    except Exception:
        payload = {}
    payload.pop(client_id, None)
    if payload:
        _atomic_write_text(cache, json.dumps(payload, indent=2), encoding='utf-8')
        os.chmod(str(cache), 0o600)
    else:
        cache.unlink(missing_ok=True)


def clear_pubmlst_auth(remove_credentials=False):
    """Forget cached OAuth tokens and optionally the saved client credentials."""
    global _oauth_credentials
    _oauth_credentials = None
    TOKENS_FILE.unlink(missing_ok=True)
    if remove_credentials:
        CONFIG_FILE.unlink(missing_ok=True)


def connect_pubmlst(client_id=None, client_secret=None, verifier=None, auth_db=None):
    set_oauth_credentials(
        client_id=client_id,
        client_secret=client_secret,
        verifier=verifier,
        auth_db=auth_db,
    )
    if not _oauth_credentials:
        raise RuntimeError('PubMLST credentials are missing. Provide client_id/client_secret.')
    client = _build_authenticated_session()
    save_pubmlst_client_credentials(
        _oauth_credentials['client_id'],
        _oauth_credentials['client_secret'],
        _oauth_credentials['auth_db'],
    )
    return client


def _exchange_client_credentials_for_access_tokens(
        client_id, client_secret, verifier=None, auth_db=None):
    import requests
    from oauthlib.oauth1 import SIGNATURE_TYPE_QUERY
    from requests_oauthlib import OAuth1

    auth_db = _validate_pubmlst_database_name(auth_db or AUTH_TEST_DB)
    request_token_url = f'{REST_BASE}/db/{auth_db}/oauth/get_request_token'
    access_token_url = f'{REST_BASE}/db/{auth_db}/oauth/get_access_token'

    def signed_token_get(url, auth_factory, label):
        for attempt in range(4):
            try:
                response = requests.get(
                    url,
                    auth=auth_factory(),
                    headers={'User-Agent': 'FastMLST'},
                    timeout=60,
                    allow_redirects=False,
                )
            except (requests.RequestException, OSError) as exc:
                if attempt >= 3:
                    raise RuntimeError(
                        f'Unable to contact the PubMLST OAuth {label} endpoint '
                        f'after 4 attempts ({type(exc).__name__}).'
                    ) from None
                time.sleep(2 ** attempt)
                continue
            if response.status_code == 429 or response.status_code >= 500:
                if attempt >= 3:
                    return response
                retry_after = (getattr(response, 'headers', {}) or {}).get(
                    'Retry-After'
                )
                try:
                    delay = float(retry_after)
                except (TypeError, ValueError):
                    delay = 2 ** attempt
                close_response = getattr(response, 'close', None)
                if callable(close_response):
                    close_response()
                time.sleep(max(0, min(delay, 120)))
                continue
            return response
        raise RuntimeError(f'Unable to contact the PubMLST OAuth {label} endpoint.')

    req_response = signed_token_get(
        request_token_url,
        lambda: OAuth1(
            client_id,
            client_secret=client_secret,
            callback_uri='oob',
            signature_type=SIGNATURE_TYPE_QUERY,
        ),
        'request-token',
    )
    if req_response.status_code != 200:
        status_code = req_response.status_code
        response_text = req_response.text
        req_response.close()
        raise RuntimeError(
            f'Unable to obtain OAuth request token from PubMLST '
            f'(HTTP {status_code}): {response_text}'
        )
    try:
        req_tokens = req_response.json()
    finally:
        req_response.close()
    request_token = req_tokens.get('oauth_token')
    request_secret = req_tokens.get('oauth_token_secret')
    if not request_token or not request_secret:
        raise RuntimeError('Unable to obtain OAuth request token from PubMLST.')

    auth_url = f'{PUBMLST_AUTHORIZE_URL}?db={auth_db}&page=authorizeClient&oauth_token={request_token}'
    if not verifier:
        print('\nAuthorize FastMLST in your browser and paste the verifier code:')
        print(auth_url)
        try:
            verifier = input('OAuth verifier: ').strip()
        except EOFError as exc:
            raise RuntimeError(
                'PubMLST authorization requires an interactive verifier. Run '
                'fastmlst --pubmlst-connect in a terminal.'
            ) from exc
    if not verifier:
        raise RuntimeError('OAuth verifier was not provided.')

    access_response = signed_token_get(
        access_token_url,
        lambda: OAuth1(
            client_id,
            client_secret=client_secret,
            resource_owner_key=request_token,
            resource_owner_secret=request_secret,
            verifier=verifier,
            signature_type=SIGNATURE_TYPE_QUERY,
        ),
        'access-token',
    )
    if access_response.status_code != 200:
        status_code = access_response.status_code
        response_text = access_response.text
        access_response.close()
        raise RuntimeError(
            f'Unable to obtain OAuth access token from PubMLST '
            f'(HTTP {status_code}): {response_text}'
        )
    try:
        access_tokens = access_response.json()
    finally:
        access_response.close()
    access_token = access_tokens.get('oauth_token')
    access_secret = access_tokens.get('oauth_token_secret')
    if not access_token or not access_secret:
        raise RuntimeError('Unable to obtain OAuth access token from PubMLST.')
    return access_token, access_secret


def _build_authenticated_session():
    if not _oauth_credentials:
        return None
    auth_db = _validate_pubmlst_database_name(
        _oauth_credentials.get('auth_db') or AUTH_TEST_DB
    )
    access_token = _oauth_credentials.get('access_token')
    access_secret = _oauth_credentials.get('access_secret')
    used_cached_token = False
    if not access_token or not access_secret:
        cached_token, cached_secret = _load_cached_access_token(
            _oauth_credentials['client_id'], auth_db
        )
        if cached_token and cached_secret:
            access_token, access_secret = cached_token, cached_secret
            used_cached_token = True
        else:
            access_token, access_secret = _exchange_client_credentials_for_access_tokens(
                _oauth_credentials['client_id'],
                _oauth_credentials['client_secret'],
                verifier=_oauth_credentials.get('verifier'),
                auth_db=auth_db,
            )
            _save_cached_access_token(
                _oauth_credentials['client_id'], access_token, access_secret, auth_db
            )

    def session_token_provider(access_session):
        token_url = f'{REST_BASE}/db/{auth_db}/oauth/get_session_token'
        response = None
        for attempt in range(4):
            try:
                response = access_session.get(
                    token_url,
                    headers={'Accept': 'application/json', 'User-Agent': 'FastMLST'},
                    timeout=60,
                    allow_redirects=False,
                )
            except OSError:
                if attempt >= 3:
                    raise
                time.sleep(2 ** attempt)
                continue
            if response.status_code == 429 or response.status_code >= 500:
                if attempt >= 3:
                    break
                retry_after = (getattr(response, 'headers', {}) or {}).get('Retry-After')
                try:
                    delay = float(retry_after)
                except (TypeError, ValueError):
                    delay = 2 ** attempt
                close_response = getattr(response, 'close', None)
                if callable(close_response):
                    close_response()
                time.sleep(max(0, min(delay, 120)))
                continue
            break
        if response.status_code in (401, 403):
            status_code = response.status_code
            close_response = getattr(response, 'close', None)
            if callable(close_response):
                close_response()
            raise PubMLSTAuthenticationError(
                f'PubMLST rejected OAuth access for {auth_db} '
                f'(HTTP {status_code}). The access token may be revoked, or the '
                'PubMLST account may not be registered for this database.'
            )
        if response.status_code >= 300:
            status_code = response.status_code
            close_response = getattr(response, 'close', None)
            if callable(close_response):
                close_response()
            raise PubMLSTHTTPError(
                f'Unable to obtain PubMLST session token (HTTP {status_code}).',
                status_code,
                token_url,
            )
        try:
            return response.json()
        except ValueError as exc:
            raise PubMLSTAuthenticationError(
                'Invalid PubMLST OAuth session token response.'
            ) from exc
        finally:
            close_response = getattr(response, 'close', None)
            if callable(close_response):
                close_response()

    def make_client(token, secret):
        client = PubMLSTClient(
            client_id=_oauth_credentials['client_id'],
            client_secret=_oauth_credentials['client_secret'],
            access_token=token,
            access_secret=secret,
            session_token_provider=session_token_provider,
            timeout=(30, 120),
            max_retries=3,
        )
        client.ensure_authenticated()
        return client

    try:
        return make_client(access_token, access_secret)
    except PubMLSTAuthenticationError:
        if not used_cached_token:
            raise
        logger.warning('Cached PubMLST access token was rejected; requesting a new one.')
        _delete_cached_access_token(_oauth_credentials['client_id'])
        access_token, access_secret = _exchange_client_credentials_for_access_tokens(
            _oauth_credentials['client_id'],
            _oauth_credentials['client_secret'],
            verifier=_oauth_credentials.get('verifier'),
            auth_db=auth_db,
        )
        _save_cached_access_token(
            _oauth_credentials['client_id'], access_token, access_secret, auth_db
        )
        return make_client(access_token, access_secret)


def fetch_json(url, session=None):
    if session is not None:
        if hasattr(session, 'get_json'):
            return session.get_json(url)
        response = session.get(url, headers={'Accept': 'application/json'})
        response.raise_for_status()
        return response.json()
    request = Request(url, headers={'Accept': 'application/json', 'User-Agent': 'FastMLST'})
    with urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode('utf-8'))


def fetch_text(url, content_type='text/plain', session=None):
    if session is not None:
        if hasattr(session, 'get_text'):
            return session.get_text(url, content_type=content_type)
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


def save_local_scheme_meta(
        codename, database, scheme_id, description, remote_snapshot,
        locus_names=None):
    p = local_scheme_meta_path(codename)
    save_local_scheme_meta_to_directory(
        p.parent,
        database,
        scheme_id,
        description,
        remote_snapshot,
        locus_names=locus_names,
    )


def save_local_scheme_meta_to_directory(
        scheme_dir, database, scheme_id, description, remote_snapshot,
        locus_names=None):
    scheme_dir = Path(scheme_dir)
    scheme_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        'database': database,
        'scheme_id': scheme_id,
        'description': description,
        'remote_snapshot': remote_snapshot,
        'loci': sorted(locus_names or []),
        'saved_at': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
    }
    _atomic_write_text(
        scheme_dir / SCHEME_LOCAL_META_FILE,
        json.dumps(payload, indent=2),
        encoding='utf-8',
    )


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
    if not profiles.is_file() or profiles.stat().st_size == 0:
        return False
    scheme_dir = profiles.parent
    for stem in stems:
        fasta = scheme_dir / f'{stem}.tfa'
        if not fasta.is_file() or fasta.stat().st_size == 0:
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
    expected_loci = set(meta.get('loci') or [])
    if expected_loci and expected_loci != set(stems):
        return False
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
    remote_locus_count = remote.get('locus_count')
    if remote_locus_count is not None:
        try:
            if len(stems) != int(remote_locus_count):
                return False
        except (TypeError, ValueError):
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
    failures = []
    unauthorized_databases = set()

    def record_failure(database, url, exc):
        failures.append(f'{url}: {exc}')
        if (
            isinstance(exc, PubMLSTAuthenticationError)
            or (
                isinstance(exc, PubMLSTHTTPError)
                and exc.status_code in (401, 403)
            )
        ):
            unauthorized_databases.add(database)

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
        try:
            schemes_payload = fetch_json(schemes_url, session=session)
        except Exception as exc:
            record_failure(db_name, schemes_url, exc)
            continue
        if not isinstance(schemes_payload, dict) or 'schemes' not in schemes_payload:
            failures.append(f'{schemes_url}: API response has no schemes field')
            continue
        schemes = schemes_payload.get('schemes')
        if not isinstance(schemes, list):
            failures.append(f'{schemes_url}: schemes field is not a list')
            continue
        if not schemes and schemes_payload.get('records') != 0:
            failures.append(f'{schemes_url}: unexpectedly empty schemes response')
            continue
        for entry in schemes:
            scheme_url = entry.get('scheme')
            if not scheme_url:
                failures.append(f'{schemes_url}: scheme entry has no scheme URL')
                continue
            try:
                detail = fetch_json(scheme_url, session=session)
            except Exception as exc:
                record_failure(db_name, scheme_url, exc)
                continue
            loci_urls = detail.get('loci') or []
            if not loci_urls:
                try:
                    loci_urls = fetch_json(
                        f'{scheme_url}/loci', session=session
                    ).get('loci', [])
                except Exception as exc:
                    record_failure(db_name, f'{scheme_url}/loci', exc)
                    continue
            description = (detail.get('description') or entry.get('description') or '')
            if not loci_urls:
                failures.append(f'{scheme_url}: API returned no loci')
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
    if failures:
        preview = '; '.join(failures[:5])
        if len(failures) > 5:
            preview += f'; ... and {len(failures) - 5} more'
        authorization_hint = ''
        if unauthorized_databases:
            database_list = ', '.join(sorted(unauthorized_databases))
            authorization_hint = (
                ' OAuth access was repeatedly rejected for: '
                f'{database_list}. This may mean that your PubMLST account is '
                'not registered for these databases. Open MY ACCOUNT > Database '
                'registrations and confirm the registrations. If they are already '
                'registered, run fastmlst --pubmlst-reset-auth, reconnect, and retry.'
            )
        raise RuntimeError(
            f'PubMLST catalog discovery was incomplete ({len(failures)} failure(s)); '
            f'the existing catalog was not replaced.{authorization_hint} {preview}'
        )
    validate_scheme_catalog(catalog)
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
        if not re.fullmatch(r'[A-Za-z0-9_]+', db) or not sid.isdigit():
            return None, None
        return db, int(sid)
    m = re.fullmatch(r'(.+)_(\d+)', s)
    if not m:
        return None, None
    db = m.group(1).strip()
    if not re.fullmatch(r'[A-Za-z0-9_]+', db):
        return None, None
    return db, int(m.group(2))


def fetch_scheme_item_direct(database, scheme_id, session=None):
    """Load one scheme from the API without scanning all databases."""
    database = _validate_pubmlst_database_name(database)
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


def validate_scheme_catalog(catalog, require_download_fields=True):
    """Reject empty, malformed, or duplicate catalogs before publication."""
    if not isinstance(catalog, list) or not catalog:
        raise ValueError('PubMLST scheme catalog is empty or invalid.')
    seen_keys = set()
    seen_codenames = set()
    required = ['codename', 'database', 'scheme_id', 'description']
    if require_download_fields:
        required.extend(('profiles_csv', 'loci'))
    for position, item in enumerate(catalog, start=1):
        if not isinstance(item, dict):
            raise ValueError(f'Catalog entry {position} is not an object.')
        missing = [
            key for key in required
            if key not in item or (key != 'description' and not item.get(key))
        ]
        if missing:
            raise ValueError(
                f'Catalog entry {position} is missing: {", ".join(missing)}.'
            )
        key = (item['database'], int(item['scheme_id']))
        if key in seen_keys or item['codename'] in seen_codenames:
            raise ValueError(
                f'Duplicate scheme in PubMLST catalog: {item["codename"]}.'
            )
        seen_keys.add(key)
        seen_codenames.add(item['codename'])
    return catalog


def save_scheme_catalog(catalog):
    validate_scheme_catalog(catalog)
    pathdb.mkdir(parents=True, exist_ok=True)
    _atomic_write_text(
        pathdb / SCHEME_CATALOG_FILE,
        json.dumps(catalog, indent=2),
        encoding='utf-8',
    )
    meta = {
        'updated_at': datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        'count': len(catalog),
    }
    _atomic_write_text(
        pathdb / SCHEME_CATALOG_META_FILE,
        json.dumps(meta, indent=2),
        encoding='utf-8',
    )


def load_scheme_catalog():
    catalog_path = pathdb / SCHEME_CATALOG_FILE
    if catalog_path.is_file():
        try:
            catalog = json.loads(catalog_path.read_text(encoding='utf-8'))
            return validate_scheme_catalog(catalog, require_download_fields=False)
        except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
            logger.warning('Ignoring invalid cached scheme catalog %s: %s', catalog_path, exc)
    bundled = _package_bundle_dir() / SCHEME_CATALOG_FILE
    if bundled.is_file():
        catalog = json.loads(bundled.read_text(encoding='utf-8'))
        return validate_scheme_catalog(catalog, require_download_fields=False)
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
    if force_refresh:
        catalog = fetch_public_scheme_catalog(session=session)
        with _database_update_lock(pathdb):
            _recover_interrupted_database_swap(pathdb)
            save_scheme_catalog(catalog)
        return catalog, False, True
    cached = load_scheme_catalog()
    if cached:
        return cached, True, False
    catalog = fetch_public_scheme_catalog(session=session)
    with _database_update_lock(pathdb):
        _recover_interrupted_database_swap(pathdb)
        save_scheme_catalog(catalog)
    return catalog, False, True


def print_remote_scheme_catalog(
        catalog, from_cache=False, meta=None, refresh_command=None):
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
    print(f'Total remote schemes: {len(catalog)}{suffix}')
    if from_cache:
        command = refresh_command or 'fastmlst --scheme-list-update'
        print(
            'Refresh the available-scheme catalog from PubMLST '
            f'(OAuth required):\n  {command}'
        )
    print()
    for i, item in enumerate(catalog, start=1):
        print(
            f"({i}) {item['codename']} | {item['database']}:{item['scheme_id']} | "
            f"{item.get('species', '')} | {item['description']}"
        )


def _validate_locus_name(locus_name):
    locus_name = str(locus_name or '').strip()
    if (
        not locus_name
        or locus_name in ('.', '..')
        or Path(locus_name).name != locus_name
        or '/' in locus_name
        or '\\' in locus_name
    ):
        raise RuntimeError(f'Unsafe or empty locus name returned by PubMLST: {locus_name!r}')
    return locus_name


def _validate_fasta_text(fasta_data, locus_name):
    if not fasta_data or not fasta_data.strip():
        raise RuntimeError(f'Empty allele FASTA returned for locus {locus_name!r}.')
    try:
        first = next(SeqIO.parse(StringIO(fasta_data), 'fasta'), None)
    except Exception as exc:
        raise RuntimeError(
            f'Invalid allele FASTA returned for locus {locus_name!r}: {exc}'
        ) from exc
    if first is None or not str(first.seq):
        raise RuntimeError(f'Allele FASTA for locus {locus_name!r} has no sequences.')


def validate_scheme_directory(scheme_dir, item, locus_names):
    """Validate a complete staged scheme before it can replace an installed copy."""
    scheme_dir = Path(scheme_dir)
    codename = item['codename']
    profile = scheme_dir / f'{codename}.txt'
    if not profile.is_file() or profile.stat().st_size == 0:
        raise RuntimeError(f'Profile table is empty for scheme {codename}.')
    profile_lines = [
        line for line in profile.read_text(encoding='utf-8', errors='replace').splitlines()
        if line.strip() and not line.lstrip().startswith('#')
    ]
    if not profile_lines or '\t' not in profile_lines[0]:
        raise RuntimeError(f'Profile table header is invalid for scheme {codename}.')

    expected_count = len(item.get('loci') or [])
    remote_count = (item.get('remote_version') or {}).get('locus_count')
    if remote_count is not None:
        try:
            remote_count = int(remote_count)
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f'Invalid locus_count for scheme {codename}.') from exc
        if remote_count != expected_count:
            raise RuntimeError(
                f'Catalog mismatch for scheme {codename}: {expected_count} locus URLs '
                f'but locus_count={remote_count}.'
            )

    expected_names = set(locus_names)
    if len(locus_names) != expected_count or len(expected_names) != expected_count:
        raise RuntimeError(
            f'Incomplete or duplicate locus download for scheme {codename}: '
            f'{len(expected_names)} of {expected_count}.'
        )
    profile_columns = [column.strip() for column in profile_lines[0].split('\t')]
    missing_profile_loci = sorted(expected_names - set(profile_columns))
    if missing_profile_loci:
        raise RuntimeError(
            f'Profile table for {codename} is missing locus column(s): '
            f'{missing_profile_loci[:10]}.'
        )
    if len(profile_columns) != len(set(profile_columns)):
        raise RuntimeError(f'Profile table for {codename} has duplicate columns.')
    installed_names = {p.stem for p in scheme_dir.glob('*.tfa')}
    if installed_names != expected_names:
        missing = sorted(expected_names - installed_names)
        unexpected = sorted(installed_names - expected_names)
        raise RuntimeError(
            f'Locus files do not match scheme {codename}; missing={missing[:5]}, '
            f'unexpected={unexpected[:5]}.'
        )
    for locus_name in sorted(expected_names):
        fasta = scheme_dir / f'{locus_name}.tfa'
        if not fasta.is_file() or fasta.stat().st_size == 0:
            raise RuntimeError(f'Empty allele file for {codename}/{locus_name}.')
        with fasta.open('r', encoding='utf-8', errors='replace') as handle:
            if next(SeqIO.parse(handle, 'fasta'), None) is None:
                raise RuntimeError(f'Invalid allele file for {codename}/{locus_name}.')
    return True


def _replace_directory_with_rollback(staged, destination):
    staged = Path(staged)
    destination = Path(destination)
    backup = destination.parent / f'.{destination.name}.backup-{uuid.uuid4().hex}'
    moved_old = False
    try:
        if destination.exists():
            os.replace(str(destination), str(backup))
            moved_old = True
        os.replace(str(staged), str(destination))
    except Exception:
        if moved_old and backup.exists() and not destination.exists():
            os.replace(str(backup), str(destination))
        raise
    if backup.exists():
        shutil.rmtree(str(backup), ignore_errors=True)


def download_scheme_data(items):
    """
    Download one scheme. Never raises: returns
    (ok, codename, loci_names, description, failed_item, error_message).
    Downloads to a temporary directory and preserves the installed copy on failure.
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
    staged_outdir = None
    try:
        schemes_root = pathdb / 'schemes'
        schemes_root.mkdir(exist_ok=True, parents=True)
        staged_outdir = Path(tempfile.mkdtemp(
            prefix=f'.{codename}.download-', dir=str(schemes_root)
        ))
        expected_files = 1 + len(item.get('loci') or [])
        bytes_downloaded = 0
        progress = tqdm(
            total=expected_files,
            desc=f'Downloading files for {codename}',
            unit='file',
            leave=False,
        )

        version_before = fetch_scheme_remote_version(
            item['database'], item['scheme_id'], session=session
        )
        before_locus_count = version_before.get('locus_count')
        if before_locus_count is not None and int(before_locus_count) != len(item['loci']):
            raise RuntimeError(
                f'Scheme {label} changed after the catalog snapshot '
                f'({len(item["loci"])} vs {before_locus_count} loci). '
                'Run --scheme-list-update and retry.'
            )

        profiles = fetch_text(
            item['profiles_csv'],
            content_type='text/tab-separated-values',
            session=session,
        )
        (staged_outdir / f'{codename}.txt').write_text(profiles, encoding='utf-8')
        bytes_downloaded += _text_nbytes(profiles)
        progress.update(1)
        progress.set_postfix_str(f'{bytes_downloaded / (1024 * 1024):.2f} MiB')

        loci_names = []
        loci_without_fasta = 0
        for locus_url in item['loci']:
            locus_record = fetch_json(locus_url, session=session)
            locus_name = _validate_locus_name(
                locus_record.get('id') or infer_locus_name_from_url(locus_url)
            )
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
            _validate_fasta_text(fasta_data, locus_name)
            (staged_outdir / f'{locus_name}.tfa').write_text(
                fasta_data, encoding='utf-8'
            )
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

        snap = fetch_scheme_remote_version(
            item['database'], item['scheme_id'], session=session
        )
        current_locus_count = snap.get('locus_count')
        if current_locus_count is not None and int(current_locus_count) != len(item['loci']):
            raise RuntimeError(
                f'Scheme {label} changed after the catalog snapshot '
                f'({len(item["loci"])} vs {current_locus_count} loci). '
                'Run --scheme-list-update and retry.'
            )
        version_keys = ('records', 'locus_count', 'last_updated', 'last_added')
        if any(
            version_before.get(key) != snap.get(key)
            for key in version_keys
            if version_before.get(key) is not None or snap.get(key) is not None
        ):
            raise RuntimeError(
                f'Scheme {label} changed while it was being downloaded; retry the update.'
            )
        validate_scheme_directory(staged_outdir, item, loci_names)
        save_local_scheme_meta_to_directory(
            staged_outdir,
            item['database'],
            item['scheme_id'],
            item['description'],
            snap,
            locus_names=loci_names,
        )
        _replace_directory_with_rollback(staged_outdir, outdir)
        staged_outdir = None

        return (True, codename, loci_names, item['description'], None, None)
    except Exception as exc:
        msg = str(exc)
        logger.error('Scheme download failed for %s (%s): %s', label, codename, msg)
        if staged_outdir is not None and staged_outdir.is_dir():
            shutil.rmtree(staged_outdir, ignore_errors=True)
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
    if present; otherwise fetch from the API. Publication happens later inside
    the database transaction, never in the active tree.
    Returns (catalog, used_cache).
    """
    cached = load_scheme_catalog()
    if cached:
        try:
            validate_scheme_catalog(cached, require_download_fields=True)
            return cached, True
        except ValueError as exc:
            logger.warning('Cached catalog cannot drive downloads; refreshing it: %s', exc)
    catalog = fetch_public_scheme_catalog(session=session)
    return catalog, False


def _is_large_scheme(item):
    description = item.get('description', '')
    if scheme_type_from_description(description) in ('cgMLST', 'wgMLST'):
        return True
    value = (item.get('remote_version') or {}).get('locus_count')
    if value is None:
        value = len(item.get('loci') or [])
    try:
        return int(value) > LARGE_SCHEME_LOCUS_THRESHOLD
    except (TypeError, ValueError):
        return len(item.get('loci') or []) > LARGE_SCHEME_LOCUS_THRESHOLD


def _filter_full_catalog(catalog, include_large_schemes=False):
    if include_large_schemes:
        return list(catalog), []
    selected = [item for item in catalog if not _is_large_scheme(item)]
    excluded = [item for item in catalog if _is_large_scheme(item)]
    if excluded:
        logger.warning(
            'Skipping %d cgMLST/wgMLST or >%d-locus scheme(s) during ALL. '
            'Use --include-large-schemes to download them explicitly.',
            len(excluded),
            LARGE_SCHEME_LOCUS_THRESHOLD,
        )
    return selected, excluded


def _load_installed_scheme_state(preserve_existing):
    scheme_number = defaultdict(list)
    scheme_list = {}
    if preserve_existing and (pathdb / 'scheme_number.pkl').is_file():
        try:
            scheme_number.update(load_obj(str(pathdb / 'scheme_number.pkl')))
        except Exception as exc:
            logger.warning('Reconstructing invalid scheme_number.pkl: %s', exc)
    if preserve_existing and (pathdb / SCHEME_LIST_FILE).is_file():
        try:
            scheme_list.update(load_obj(str(pathdb / SCHEME_LIST_FILE)))
        except Exception as exc:
            logger.warning('Reconstructing invalid scheme_list.pkl: %s', exc)
    schemes_root = pathdb / 'schemes'
    if preserve_existing and schemes_root.is_dir():
        for scheme_dir in sorted(p for p in schemes_root.iterdir() if p.is_dir()):
            if scheme_dir.name.startswith('.'):
                continue
            codename = scheme_dir.name
            loci = sorted(p.stem for p in scheme_dir.glob('*.tfa'))
            profile = scheme_dir / f'{codename}.txt'
            if loci and profile.is_file():
                if codename not in scheme_number:
                    meta = load_local_scheme_meta(codename) or {}
                    expected = sorted(set(meta.get('loci') or []))
                    scheme_number[codename] = expected or loci
                meta = load_local_scheme_meta(codename) or {}
                scheme_list.setdefault(codename, meta.get('description', ''))
    return scheme_number, scheme_list


def _validate_and_clean_installed_manifest(scheme_number):
    """Validate manifest inputs and remove stale loci inside the staging tree."""
    schemes_root = pathdb / 'schemes'
    for codename, loci in scheme_number.items():
        scheme_dir = schemes_root / codename
        profile = scheme_dir / f'{codename}.txt'
        if not profile.is_file() or profile.stat().st_size == 0:
            raise RuntimeError(f'Missing profile table required by manifest: {profile}')
        expected = set(loci)
        if not expected:
            raise RuntimeError(f'Scheme {codename} has no loci in the manifest.')
        profile_header = next(
            (
                line for line in profile.read_text(
                    encoding='utf-8', errors='replace'
                ).splitlines()
                if line.strip() and not line.lstrip().startswith('#')
            ),
            '',
        )
        profile_columns = {column.strip() for column in profile_header.split('\t')}
        missing_profile_loci = sorted(expected - profile_columns)
        if missing_profile_loci:
            raise RuntimeError(
                f'Profile table for {codename} is missing manifest loci: '
                f'{missing_profile_loci[:10]}'
            )
        installed = {p.stem: p for p in scheme_dir.glob('*.tfa')}
        missing = sorted(expected - set(installed))
        if missing:
            raise RuntimeError(
                f'Scheme {codename} is incomplete; missing loci: {missing[:10]}'
            )
        for stale in sorted(set(installed) - expected):
            logger.warning('Removing stale locus from staged scheme: %s/%s', codename, stale)
            installed[stale].unlink()


def _iter_database_fasta_records(database_root, scheme_number):
    schemes_root = Path(database_root) / 'schemes'
    for scheme in sorted(scheme_number):
        if scheme_excluded_from_mlst_concat(scheme):
            logger.info('Omitted from mlst.fasta merge: %s', scheme)
            continue
        for locus in sorted(set(scheme_number[scheme])):
            fasta = schemes_root / scheme / f'{locus}.tfa'
            if not fasta.is_file() or fasta.stat().st_size == 0:
                raise RuntimeError(f'Missing allele file required by manifest: {fasta}')
            parsed_any = False
            for record in SeqIO.parse(str(fasta), 'fasta'):
                parsed_any = True
                record.id = f'{scheme}.{record.id}'
                record.description = ''
                yield record
            if not parsed_any:
                raise RuntimeError(f'Invalid or empty allele FASTA: {fasta}')


def _blast_index_paths(database_root, prefix='mlst.fasta'):
    database_root = Path(database_root)
    return sorted(
        p for p in database_root.glob(f'{prefix}.*')
        if p.is_file() and any(p.name.endswith(suffix) for suffix in BLAST_INDEX_SUFFIXES)
    )


def blast_database_is_ready(database_root=None, deep=False):
    database_root = Path(database_root or pathdb)
    fasta = database_root / 'mlst.fasta'
    if not fasta.is_file() or fasta.stat().st_size == 0:
        return False
    indexes = _blast_index_paths(database_root)
    names = [p.name for p in indexes]
    ready = all(
        any(name.endswith(suffix) for name in names)
        for suffix in ('.nhr', '.nin', '.nsq')
    )
    if not ready or not deep:
        return ready
    try:
        result = subprocess.run(
            ['blastdbcmd', '-db', str(database_root / 'mlst.fasta'), '-info'],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.returncode == 0


def _build_database_artifacts(
        scheme_number, scheme_list,
        resource_source='https://rest.pubmlst.org/db',
        blast_title_prefix='PubMLST'):
    if not scheme_number:
        raise RuntimeError('No MLST schemes were available after update.')
    pathdb.mkdir(parents=True, exist_ok=True)
    _validate_and_clean_installed_manifest(scheme_number)
    _atomic_write_text(
        pathdb / DBASES_COMPAT_FILE,
        f'<resources source="{resource_source}"/>',
        encoding='utf-8',
    )
    save_obj(dict(scheme_number), pathdb / 'scheme_number.pkl')
    save_obj(dict(scheme_list), pathdb / SCHEME_LIST_FILE)

    build_id = uuid.uuid4().hex
    temporary_fasta = pathdb / f'.mlst.fasta-{build_id}.tmp'
    final_prefix = pathdb / 'mlst.fasta'
    try:
        count = SeqIO.write(
            _iter_database_fasta_records(pathdb, scheme_number),
            str(temporary_fasta),
            'fasta',
        )
        if count == 0:
            raise RuntimeError(
                'No standard MLST allele sequences were available for the BLAST database. '
                'Install at least one non-rMLST scheme.'
            )
        # The BLAST v5 .njs metadata embeds the database prefix. Build with the
        # final prefix inside the isolated staging tree; renaming these files
        # afterwards would leave an unusable database.
        for old_index in _blast_index_paths(pathdb):
            old_index.unlink()
        blastdb_cmd = [
            'makeblastdb',
            '-hash_index',
            '-in', str(temporary_fasta),
            '-out', str(final_prefix),
            '-dbtype', 'nucl',
            '-title', f'{blast_title_prefix}_{date.today().strftime("%d%m%y")}',
            '-parse_seqids',
        ]
        result = subprocess.run(
            blastdb_cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f'makeblastdb failed while building {blast_title_prefix} database: '
                + (result.stderr or '').strip()
            )
        if not all(
            any(p.name.endswith(suffix) for p in _blast_index_paths(pathdb))
            for suffix in ('.nhr', '.nin', '.nsq')
        ):
            raise RuntimeError('makeblastdb returned success without complete index files.')

        os.replace(str(temporary_fasta), str(final_prefix))
        if not blast_database_is_ready(pathdb):
            raise RuntimeError('Published BLAST database failed validation.')
        validation = subprocess.run(
            ['blastdbcmd', '-db', str(final_prefix), '-info'],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if validation.returncode != 0:
            raise RuntimeError(
                f'blastdbcmd could not open the generated {blast_title_prefix} database: '
                + (validation.stderr or '').strip()
            )
    finally:
        temporary_fasta.unlink(missing_ok=True)


def _update_database_in_place(
        scheme_catalog, session, threads, preserve_existing,
        force_redownload_schemes=False):
    scheme_number, scheme_list = _load_installed_scheme_state(preserve_existing)
    to_download = []
    if force_redownload_schemes:
        to_download = list(scheme_catalog)
    else:
        for item in scheme_catalog:
            codename = item['codename']
            if local_scheme_matches_remote(
                codename, item['database'], item['scheme_id'], session=session
            ):
                logger.info('Skipping unchanged scheme %s', codename)
                scheme_number[codename] = local_scheme_locus_stems(codename)
                meta = load_local_scheme_meta(codename) or {}
                scheme_list[codename] = meta.get('description') or item.get('description', '')
            else:
                to_download.append(item)

    failures = []
    if to_download:
        with ThreadPool(threads) as pool:
            tasks = [(item, session) for item in to_download]
            for result in tqdm(
                pool.imap(download_scheme_data, tasks),
                total=len(tasks),
                desc=f'Downloading schemes using {threads} thread(s)',
                unit='scheme',
                leave=True,
            ):
                ok, codename, loci_names, description, failed_item, error = result
                if ok:
                    scheme_number[codename] = loci_names
                    scheme_list[codename] = description
                else:
                    failures.append((failed_item, error))
    if failures:
        details = '; '.join(
            f"{item.get('database')}:{item.get('scheme_id')}: {error}"
            for item, error in failures[:10]
        )
        raise RuntimeError(
            f'{len(failures)} requested PubMLST scheme download(s) failed; '
            f'no changes were published. {details}'
        )
    _build_database_artifacts(scheme_number, scheme_list)


def update_mlstdb_selected(
        threads, selectors=None, force_redownload_schemes=False,
        include_large_schemes=False):
    """Build a verified database tree and publish it as one transaction."""
    try:
        threads = int(threads)
    except (TypeError, ValueError) as exc:
        raise ValueError('threads must be an integer greater than zero.') from exc
    if threads < 1:
        raise ValueError('threads must be an integer greater than zero.')
    active_database = _validate_database_destination(pathdb)
    if not _oauth_credentials:
        configure_pubmlst_auth_from_environment(require_auth=True)
    with _database_update_lock(active_database):
        _recover_interrupted_database_swap(active_database)
        session = _build_authenticated_session()
        if session is None:
            raise RuntimeError('Authenticated PubMLST access is required for downloads.')
        logger.info('Using authenticated PubMLST API access (OAuth1)')
        if threads > 1:
            logger.info('Using one download thread for the shared OAuth session.')
            threads = 1
        stage = None
        try:
            if selectors:
                catalog = build_scheme_catalog_direct(selectors, session=session)
                download_catalog = catalog
                preserve_existing = True
            else:
                logger.info('Refreshing the complete PubMLST catalog before ALL update.')
                catalog = fetch_public_scheme_catalog(session=session)
                download_catalog, _excluded = _filter_full_catalog(
                    catalog, include_large_schemes=include_large_schemes
                )
                preserve_existing = False
            if not download_catalog:
                raise RuntimeError('No PubMLST schemes matched the requested update.')

            stage = _prepare_database_staging(
                active_database, clone_existing=preserve_existing
            )
            with _using_database_path(stage):
                if not preserve_existing:
                    save_scheme_catalog(catalog)
                _update_database_in_place(
                    download_catalog,
                    session,
                    threads,
                    preserve_existing=preserve_existing,
                    force_redownload_schemes=force_redownload_schemes,
                )
            _publish_database_staging(active_database, stage)
            stage = None
            logger.info('PubMLST update completed and published atomically.')
        finally:
            if stage is not None and stage.exists():
                shutil.rmtree(str(stage), ignore_errors=True)
            close = getattr(session, 'close', None)
            if callable(close):
                close()

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
