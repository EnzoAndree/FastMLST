import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from fastmlst import update_mlst_kit as update


BLAST_AVAILABLE = all(
    shutil.which(command) for command in ('makeblastdb', 'blastdbcmd', 'blastn')
)
FIXTURE_ROOT = Path(__file__).parent / 'fixtures' / 'cdifficile_r20291'
GENOME = FIXTURE_ROOT / 'Clostridioides_difficile_R20291_FN545816.1.fasta.gz'
SCHEME = 'pubmlst_cdifficile_seqdef_1'
LOCI = ['adk', 'atpA', 'dxr', 'glyA', 'recA', 'sodA', 'tpi']
EXPECTED_RESULT = (
    'Clostridioides_difficile_R20291_FN545816.1.fasta.gz,'
    'pubmlst_cdifficile_seqdef_1,1,adk(1),atpA(1),dxr(1),glyA(10),'
    'recA(1),sodA(3),tpi(5),mlst_clade(2)'
)


@unittest.skipUnless(BLAST_AVAILABLE, 'NCBI BLAST+ executables are not installed')
class ClostridioidesDifficileR20291IntegrationTests(unittest.TestCase):
    def test_cli_assigns_expected_scheme_st_and_alleles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            database = Path(tmpdir) / 'pubmlst'
            destination = database / 'schemes' / SCHEME
            shutil.copytree(FIXTURE_ROOT / 'schemes' / SCHEME, destination)

            with update._using_database_path(database):
                update._build_database_artifacts(
                    {SCHEME: LOCI},
                    {SCHEME: 'MLST'},
                )

            result = subprocess.run(
                [
                    sys.executable,
                    '-m',
                    'fastmlst.cli',
                    '--db_path',
                    str(database),
                    '--threads',
                    '1',
                    str(GENOME),
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=120,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout.strip(), EXPECTED_RESULT)


if __name__ == '__main__':
    unittest.main()
