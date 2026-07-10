import gzip
import shutil
import subprocess
import tempfile
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from Bio import SeqIO

from fastmlst import cli
from fastmlst import ngstar
from fastmlst import update_mlst_kit as database_kit


FIXTURE_ROOT = Path(__file__).parent / 'fixtures' / 'ngstar_v2_source'
FASTMLST_EXECUTABLE = shutil.which('fastmlst')
BLAST_AVAILABLE = all(
    shutil.which(command) for command in ('makeblastdb', 'blastdbcmd', 'blastn')
) and FASTMLST_EXECUTABLE is not None
EXPECTED_LOCUS_COUNTS = {
    'penA': 962,
    'mtrR': 889,
    'porB': 87,
    'ponA': 23,
    'gyrA': 84,
    'parC': 243,
    '23S': 112,
}
ST1_ALLELES = {
    'penA': '22.001',
    'mtrR': '10.0',
    'porB': '13.0',
    'ponA': '100.0',
    'gyrA': '100.0',
    'parC': '2.0',
    '23S': '100.0',
}


class FakeResponse:
    def __init__(self, content, status_code=200, headers=None):
        self.content = content
        self.status_code = status_code
        self.headers = headers or {}
        self.closed = False

    def close(self):
        self.closed = True


class FixtureSession:
    def __init__(self, invalid_profiles=False):
        self.invalid_profiles = invalid_profiles
        self.calls = []

    def get(self, url, **kwargs):
        self.calls.append((url, kwargs))
        if url == ngstar.NGSTAR_PROFILE_URL:
            if self.invalid_profiles:
                return FakeResponse(b'not an xlsx workbook')
            return FakeResponse((FIXTURE_ROOT / 'ngstar_profiles.xlsx').read_bytes())
        for locus in ngstar.NGSTAR_LOCI:
            if url == ngstar.NGSTAR_ALLELE_URL.format(locus=locus):
                return FakeResponse(
                    gzip.decompress(
                        (FIXTURE_ROOT / 'alleles' / f'{locus}.fasta.gz').read_bytes()
                    )
                )
        raise AssertionError(f'Unexpected NG-STAR URL: {url}')


def source_allele_identifiers():
    identifiers = {}
    counts = {}
    for locus in ngstar.NGSTAR_LOCI:
        payload = gzip.decompress(
            (FIXTURE_ROOT / 'alleles' / f'{locus}.fasta.gz').read_bytes()
        )
        _text, identifiers[locus], counts[locus] = ngstar._validate_allele_fasta(
            locus, payload
        )
    return identifiers, counts


def write_st1_query(destination):
    records = []
    for locus, allele in ST1_ALLELES.items():
        payload = gzip.decompress(
            (FIXTURE_ROOT / 'alleles' / f'{locus}.fasta.gz').read_bytes()
        ).decode('utf-8')
        expected_id = f'{locus}_{allele}'
        record = next(
            (
                candidate
                for candidate in SeqIO.parse(StringIO(payload), 'fasta')
                if candidate.id == expected_id
            ),
            None,
        )
        if record is None:
            raise AssertionError(f'Missing fixture allele {expected_id}')
        record.id = f'contig_{locus}'
        record.name = record.id
        record.description = ''
        records.append(record)
    SeqIO.write(records, str(destination), 'fasta')


class NGStarV2SourceTests(unittest.TestCase):
    def test_cli_update_flag_uses_separate_ngstar_database_path(self):
        old_pathdb = database_kit.pathdb
        with tempfile.TemporaryDirectory() as tmpdir:
            database = Path(tmpdir) / 'NG-STAR-v2'
            argv = [
                'fastmlst',
                '--update-ngstar-v2',
                '--ngstar-db-path',
                str(database),
            ]
            try:
                with patch('sys.argv', argv):
                    with patch.object(
                        ngstar,
                        'update_ngstar_v2_database',
                        return_value={'profile_count': 8125},
                    ) as update:
                        with patch('sys.stdout', StringIO()) as output:
                            with self.assertRaises(SystemExit) as caught:
                                cli.main()
                self.assertEqual(caught.exception.code, None)
                update.assert_called_once_with(database)
                self.assertIn('8125 profiles', output.getvalue())
                self.assertEqual(database_kit.pathdb, database)
            finally:
                database_kit.pathdb = old_pathdb

    def test_official_xlsx_converts_to_exact_fasta_allele_identifiers(self):
        identifiers, counts = source_allele_identifiers()
        profiles, profile_count = ngstar.convert_ngstar_profiles(
            (FIXTURE_ROOT / 'ngstar_profiles.xlsx').read_bytes(),
            identifiers,
        )

        lines = profiles.splitlines()
        self.assertEqual(counts, EXPECTED_LOCUS_COUNTS)
        self.assertEqual(profile_count, 8125)
        self.assertEqual(lines[0], 'ST\tpenA\tmtrR\tporB\tponA\tgyrA\tparC\t23S')
        self.assertEqual(
            lines[1],
            '1\t22.001\t10.0\t13.0\t100.0\t100.0\t2.0\t100.0',
        )
        self.assertEqual(
            lines[-1],
            '8306\t5.002\t597.0\t1.0\t1.0\t7.0\t3.0\t2.0',
        )

    def test_failed_profile_conversion_preserves_installed_database(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            database = Path(tmpdir) / 'NG-STAR-v2'
            database.mkdir()
            sentinel = database / 'sentinel.txt'
            sentinel.write_text('active database', encoding='utf-8')

            with self.assertRaisesRegex(RuntimeError, 'invalid XLSX'):
                ngstar.update_ngstar_v2_database(
                    database, session=FixtureSession(invalid_profiles=True)
                )

            self.assertEqual(sentinel.read_text(encoding='utf-8'), 'active database')
            self.assertFalse(
                list(database.parent.glob(f'.{database.name}.staging-*'))
            )


@unittest.skipUnless(BLAST_AVAILABLE, 'NCBI BLAST+ executables are not installed')
class NGStarV2DatabaseIntegrationTests(unittest.TestCase):
    def test_official_snapshot_builds_independent_queryable_database(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            database = Path(tmpdir) / 'NG-STAR-v2'
            session = FixtureSession()

            metadata = ngstar.update_ngstar_v2_database(
                database, session=session
            )

            self.assertEqual(metadata['profile_count'], 8125)
            self.assertEqual(metadata['locus_counts'], EXPECTED_LOCUS_COUNTS)
            self.assertEqual(len(session.calls), 8)
            self.assertTrue(database_kit.blast_database_is_ready(database, deep=True))
            self.assertTrue((database / 'source' / 'ngstar_profiles.xlsx').is_file())
            self.assertFalse((database / 'scheme_catalog.json').exists())
            self.assertIn(
                'https://ngstar.canada.ca',
                (database / 'dbases.xml').read_text(encoding='utf-8'),
            )
            profile = (
                database
                / 'schemes'
                / ngstar.NGSTAR_SCHEME
                / f'{ngstar.NGSTAR_SCHEME}.txt'
            )
            self.assertEqual(
                profile.read_text(encoding='utf-8').splitlines()[1],
                '1\t22.001\t10.0\t13.0\t100.0\t100.0\t2.0\t100.0',
            )

            query = Path(tmpdir) / 'ngstar_st1.fasta'
            write_st1_query(query)
            result = subprocess.run(
                [
                    FASTMLST_EXECUTABLE,
                    '--ngstar-v2',
                    '--ngstar-db-path',
                    str(database),
                    '--threads',
                    '1',
                    str(query),
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=120,
                cwd=str(Path(__file__).resolve().parents[1]),
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                result.stdout.strip(),
                'ngstar_st1.fasta,ngstar_v2,1,penA(22.001),mtrR(10.0),'
                'porB(13.0),ponA(100.0),gyrA(100.0),parC(2.0),23S(100.0)',
            )


if __name__ == '__main__':
    unittest.main()
