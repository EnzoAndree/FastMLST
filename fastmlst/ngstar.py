"""Download and build an independent NG-STAR v2 typing database."""

import hashlib
import json
import logging
import os
import re
import shutil
import tempfile
import time
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from io import BytesIO, StringIO
from pathlib import Path
from urllib.parse import urlparse
from xml.etree import ElementTree
from zipfile import BadZipFile, ZipFile

import requests
from Bio import SeqIO
from tqdm import tqdm

from fastmlst import update_mlst_kit as database_kit


logger = logging.getLogger('ngstar')
NGSTAR_BASE = 'https://ngstar.canada.ca'
NGSTAR_HOST = 'ngstar.canada.ca'
NGSTAR_SCHEME = 'ngstar_v2'
NGSTAR_DESCRIPTION = 'NG-STAR v2 antimicrobial resistance sequence typing'
NGSTAR_LOCI = ('penA', 'mtrR', 'porB', 'ponA', 'gyrA', 'parC', '23S')
NGSTAR_PROFILE_URL = f'{NGSTAR_BASE}/sequence_types/download?lang=en'
NGSTAR_ALLELE_URL = (
    f'{NGSTAR_BASE}/alleles/download?lang=en&loci_name={{locus}}'
)
NGSTAR_META_FILE = 'ngstar_v2_meta.json'
MAX_XLSX_UNCOMPRESSED_BYTES = 50 * 1024 * 1024
_SPREADSHEET_NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
_DOCUMENT_REL_NS = (
    'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
)
_PACKAGE_REL_NS = (
    'http://schemas.openxmlformats.org/package/2006/relationships'
)


def default_ngstar_database_path():
    return Path.home() / '.cache' / 'fastmlst' / 'NG-STAR-v2'


def _atomic_write_bytes(destination, payload):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f'.{destination.name}.',
        suffix='.tmp',
        dir=str(destination.parent),
    )
    try:
        with os.fdopen(fd, 'wb') as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, str(destination))
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def _download_bytes(url, session, max_retries=3, sleep=time.sleep):
    parsed = urlparse(url)
    if parsed.scheme != 'https' or parsed.hostname != NGSTAR_HOST:
        raise ValueError(f'Refusing to download NG-STAR data from {url!r}.')

    for attempt in range(max_retries + 1):
        try:
            response = session.get(
                url,
                headers={
                    'Accept': '*/*',
                    'User-Agent': 'FastMLST',
                },
                timeout=(30, 180),
                allow_redirects=False,
            )
        except (requests.RequestException, OSError) as exc:
            if attempt >= max_retries:
                raise RuntimeError(
                    f'Unable to download NG-STAR data after {max_retries} '
                    f'retries ({type(exc).__name__}).'
                ) from None
            sleep(2 ** attempt)
            continue

        status_code = int(response.status_code)
        if status_code == 429 or status_code >= 500:
            if attempt >= max_retries:
                response.close()
                raise RuntimeError(
                    f'NG-STAR download failed after {max_retries} retries '
                    f'(HTTP {status_code}).'
                )
            retry_after = (getattr(response, 'headers', {}) or {}).get(
                'Retry-After'
            )
            try:
                delay = float(retry_after)
            except (TypeError, ValueError):
                delay = 2 ** attempt
            response.close()
            sleep(max(0, min(delay, 120)))
            continue
        if status_code != 200:
            response.close()
            raise RuntimeError(
                f'NG-STAR download failed with HTTP {status_code}.'
            )
        try:
            payload = bytes(response.content)
        finally:
            response.close()
        if not payload:
            raise RuntimeError(f'NG-STAR returned an empty response for {url}.')
        return payload

    raise RuntimeError(f'Unable to download NG-STAR data from {url}.')


def _validate_allele_fasta(locus, payload):
    try:
        text = payload.decode('utf-8')
    except UnicodeDecodeError as exc:
        raise RuntimeError(f'NG-STAR returned invalid UTF-8 FASTA for {locus}.') from exc

    records = list(SeqIO.parse(StringIO(text), 'fasta'))
    if not records:
        raise RuntimeError(f'NG-STAR returned no FASTA records for {locus}.')

    identifiers = set()
    numeric_identifiers = {}
    prefix = f'{locus}_'
    for record in records:
        if not record.id.startswith(prefix) or not str(record.seq):
            raise RuntimeError(
                f'NG-STAR returned an invalid FASTA record for {locus}: '
                f'{record.id!r}.'
            )
        if record.id in identifiers:
            raise RuntimeError(
                f'NG-STAR returned duplicate FASTA identifier {record.id!r}.'
            )
        identifiers.add(record.id)
        allele = record.id[len(prefix):]
        try:
            numeric = Decimal(allele)
        except InvalidOperation as exc:
            raise RuntimeError(
                f'NG-STAR allele identifier is not numeric: {record.id!r}.'
            ) from exc
        if numeric in numeric_identifiers:
            raise RuntimeError(
                f'NG-STAR allele identifiers {numeric_identifiers[numeric]!r} '
                f'and {allele!r} are numerically ambiguous for {locus}.'
            )
        numeric_identifiers[numeric] = allele

    if not text.endswith('\n'):
        text += '\n'
    return text, numeric_identifiers, len(records)


def _xlsx_shared_strings(archive):
    if 'xl/sharedStrings.xml' not in archive.namelist():
        return []
    root = ElementTree.fromstring(archive.read('xl/sharedStrings.xml'))
    return [
        ''.join(
            node.text or ''
            for node in item.iter(f'{{{_SPREADSHEET_NS}}}t')
        )
        for item in root.findall(f'{{{_SPREADSHEET_NS}}}si')
    ]


def _first_worksheet_path(archive):
    workbook = ElementTree.fromstring(archive.read('xl/workbook.xml'))
    sheet = workbook.find(
        f'{{{_SPREADSHEET_NS}}}sheets/{{{_SPREADSHEET_NS}}}sheet'
    )
    if sheet is None:
        raise RuntimeError('NG-STAR profile workbook has no worksheets.')
    relationship_id = sheet.get(f'{{{_DOCUMENT_REL_NS}}}id')
    relationships = ElementTree.fromstring(
        archive.read('xl/_rels/workbook.xml.rels')
    )
    target = None
    for relationship in relationships.findall(f'{{{_PACKAGE_REL_NS}}}Relationship'):
        if relationship.get('Id') == relationship_id:
            target = relationship.get('Target')
            break
    if not target:
        raise RuntimeError('NG-STAR profile worksheet relationship is missing.')
    target = target.lstrip('/')
    if not target.startswith('xl/'):
        target = f'xl/{target}'
    return target


def _column_index(reference):
    match = re.match(r'^([A-Z]+)', str(reference or '').upper())
    if not match:
        raise RuntimeError(f'Invalid XLSX cell reference: {reference!r}.')
    value = 0
    for character in match.group(1):
        value = value * 26 + ord(character) - ord('A') + 1
    return value - 1


def _cell_text(cell, shared_strings):
    cell_type = cell.get('t')
    if cell_type == 'inlineStr':
        return ''.join(
            node.text or ''
            for node in cell.iter(f'{{{_SPREADSHEET_NS}}}t')
        )
    value_node = cell.find(f'{{{_SPREADSHEET_NS}}}v')
    if value_node is None or value_node.text is None:
        return ''
    value = value_node.text.strip()
    if cell_type == 's':
        try:
            return shared_strings[int(value)]
        except (IndexError, ValueError) as exc:
            raise RuntimeError('Invalid shared string in NG-STAR workbook.') from exc
    return value


def _xlsx_rows(payload):
    try:
        with ZipFile(BytesIO(payload)) as archive:
            if sum(item.file_size for item in archive.infolist()) > MAX_XLSX_UNCOMPRESSED_BYTES:
                raise RuntimeError('NG-STAR profile workbook is unexpectedly large.')
            shared_strings = _xlsx_shared_strings(archive)
            worksheet_path = _first_worksheet_path(archive)
            worksheet = ElementTree.fromstring(archive.read(worksheet_path))
    except (BadZipFile, KeyError, ElementTree.ParseError) as exc:
        raise RuntimeError('NG-STAR returned an invalid XLSX profile workbook.') from exc

    rows = []
    for row in worksheet.findall(
            f'.//{{{_SPREADSHEET_NS}}}sheetData/{{{_SPREADSHEET_NS}}}row'):
        values = [''] * 8
        for cell in row.findall(f'{{{_SPREADSHEET_NS}}}c'):
            index = _column_index(cell.get('r'))
            if index < len(values):
                values[index] = _cell_text(cell, shared_strings).strip()
        if any(values):
            rows.append(values)
    if not rows:
        raise RuntimeError('NG-STAR profile workbook contains no rows.')
    return rows


def _canonical_sequence_type(raw_value):
    try:
        value = Decimal(raw_value)
    except InvalidOperation as exc:
        raise RuntimeError(f'Invalid NG-STAR sequence type: {raw_value!r}.') from exc
    if not value.is_finite():
        raise RuntimeError(f'Invalid NG-STAR sequence type: {raw_value!r}.')
    integral = value.to_integral_value()
    if value != integral or integral < 1:
        raise RuntimeError(f'Invalid NG-STAR sequence type: {raw_value!r}.')
    return str(int(integral))


def _canonical_allele(locus, raw_value, numeric_identifiers):
    try:
        numeric = Decimal(raw_value)
    except InvalidOperation as exc:
        raise RuntimeError(
            f'Invalid NG-STAR allele for {locus}: {raw_value!r}.'
        ) from exc
    if not numeric.is_finite():
        raise RuntimeError(
            f'Invalid NG-STAR allele for {locus}: {raw_value!r}.'
        )
    try:
        return numeric_identifiers[numeric]
    except KeyError as exc:
        raise RuntimeError(
            f'NG-STAR profile references missing {locus} allele {raw_value!r}.'
        ) from exc


def convert_ngstar_profiles(profile_xlsx, allele_identifiers):
    rows = _xlsx_rows(profile_xlsx)
    expected_header = ['Sequence Type'] + list(NGSTAR_LOCI)
    if rows[0] != expected_header:
        raise RuntimeError(
            f'Unexpected NG-STAR profile columns: {rows[0]!r}; '
            f'expected {expected_header!r}.'
        )

    output_rows = [['ST'] + list(NGSTAR_LOCI)]
    sequence_types = set()
    for row_number, row in enumerate(rows[1:], start=2):
        if len(row) != len(expected_header) or any(value == '' for value in row):
            raise RuntimeError(
                f'Incomplete NG-STAR profile workbook row {row_number}.'
            )
        sequence_type = _canonical_sequence_type(row[0])
        if sequence_type in sequence_types:
            raise RuntimeError(
                f'Duplicate NG-STAR sequence type {sequence_type!r}.'
            )
        sequence_types.add(sequence_type)
        alleles = [
            _canonical_allele(locus, row[index], allele_identifiers[locus])
            for index, locus in enumerate(NGSTAR_LOCI, start=1)
        ]
        output_rows.append([sequence_type] + alleles)

    if not sequence_types:
        raise RuntimeError('NG-STAR profile workbook contains no profiles.')
    profile_text = ''.join('\t'.join(row) + '\n' for row in output_rows)
    return profile_text, len(sequence_types)


def update_ngstar_v2_database(database_root=None, session=None):
    """Download, validate, build, and atomically publish an NG-STAR v2 database."""
    active_database = database_kit._validate_database_destination(
        database_root or default_ngstar_database_path()
    )
    owns_session = session is None
    stage = None
    with database_kit._database_update_lock(active_database):
        database_kit._recover_interrupted_database_swap(active_database)
        if session is None:
            session = requests.Session()
        try:
            stage = Path(tempfile.mkdtemp(
                prefix=f'.{active_database.name}.staging-',
                dir=str(active_database.parent),
            ))
            scheme_directory = stage / 'schemes' / NGSTAR_SCHEME
            scheme_directory.mkdir(parents=True)
            allele_identifiers = {}
            locus_counts = {}
            source_hashes = {}
            for locus in tqdm(
                    NGSTAR_LOCI,
                    desc='Downloading NG-STAR v2 alleles',
                    unit='locus'):
                url = NGSTAR_ALLELE_URL.format(locus=locus)
                payload = _download_bytes(url, session)
                text, identifiers, record_count = _validate_allele_fasta(
                    locus, payload
                )
                database_kit._atomic_write_text(
                    scheme_directory / f'{locus}.tfa', text, encoding='utf-8'
                )
                allele_identifiers[locus] = identifiers
                locus_counts[locus] = record_count
                source_hashes[f'{locus}.fasta'] = hashlib.sha256(payload).hexdigest()

            profile_xlsx = _download_bytes(NGSTAR_PROFILE_URL, session)
            profile_text, profile_count = convert_ngstar_profiles(
                profile_xlsx, allele_identifiers
            )
            database_kit._atomic_write_text(
                scheme_directory / f'{NGSTAR_SCHEME}.txt',
                profile_text,
                encoding='utf-8',
            )
            _atomic_write_bytes(
                stage / 'source' / 'ngstar_profiles.xlsx', profile_xlsx
            )
            source_hashes['ngstar_profiles.xlsx'] = hashlib.sha256(
                profile_xlsx
            ).hexdigest()

            metadata = {
                'name': 'NG-STAR',
                'version': '2.0',
                'scheme': NGSTAR_SCHEME,
                'source': NGSTAR_BASE,
                'source_urls': {
                    'profiles': NGSTAR_PROFILE_URL,
                    'alleles': {
                        locus: NGSTAR_ALLELE_URL.format(locus=locus)
                        for locus in NGSTAR_LOCI
                    },
                },
                'downloaded_at': datetime.now(timezone.utc).replace(
                    microsecond=0
                ).isoformat(),
                'profile_count': profile_count,
                'locus_counts': locus_counts,
                'sha256': source_hashes,
            }
            database_kit._atomic_write_text(
                stage / NGSTAR_META_FILE,
                json.dumps(metadata, indent=2, sort_keys=True) + '\n',
                encoding='utf-8',
            )

            with database_kit._using_database_path(stage):
                database_kit._build_database_artifacts(
                    {NGSTAR_SCHEME: list(NGSTAR_LOCI)},
                    {NGSTAR_SCHEME: NGSTAR_DESCRIPTION},
                    resource_source=NGSTAR_BASE,
                    blast_title_prefix='NG_STAR_v2',
                )
            database_kit._publish_database_staging(
                active_database, stage
            )
            stage = None
            logger.info(
                'NG-STAR v2 database updated: %s profiles under %s',
                profile_count,
                active_database,
            )
            return metadata
        finally:
            if stage is not None and stage.exists():
                shutil.rmtree(str(stage), ignore_errors=True)
            if owns_session:
                session.close()
