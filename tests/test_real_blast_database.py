import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from fastmlst import update_mlst_kit as update


BLAST_AVAILABLE = all(
    shutil.which(command) for command in ('makeblastdb', 'blastdbcmd', 'blastn')
)


@unittest.skipUnless(BLAST_AVAILABLE, 'NCBI BLAST+ executables are not installed')
class RealBlastDatabaseTests(unittest.TestCase):
    def test_generated_v5_database_uses_final_prefix_and_is_queryable(self):
        old_pathdb = update.pathdb
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                update.pathdb = Path(tmpdir) / 'pubmlst'
                scheme_dir = update.pathdb / 'schemes' / 'demo'
                scheme_dir.mkdir(parents=True)
                (scheme_dir / 'demo.txt').write_text(
                    'ST\tadk\n1\t1\n2\t2\n', encoding='utf-8'
                )
                (scheme_dir / 'adk.tfa').write_text(
                    '>adk_1\nACGTACGTACGT\n>adk_2\nACGTTCGTACGT\n', encoding='utf-8'
                )

                update._build_database_artifacts(
                    {'demo': ['adk']}, {'demo': 'MLST'}
                )

                info = subprocess.run(
                    ['blastdbcmd', '-db', str(update.pathdb / 'mlst.fasta'), '-info'],
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                self.assertEqual(info.returncode, 0, info.stderr)

                query = update.pathdb / 'query.fasta'
                query.write_text('>query\nACGTACGTACGT\n', encoding='utf-8')
                result = subprocess.run(
                    [
                        'blastn', '-query', str(query),
                        '-db', str(update.pathdb / 'mlst.fasta'),
                        '-task', 'blastn-short',
                        '-dust', 'no',
                        '-evalue', '1000',
                        '-outfmt', '6 sseqid',
                    ],
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('demo.adk_1', result.stdout)

                active_root = update.pathdb
                active_fasta = (active_root / 'mlst.fasta').read_bytes()
                stage = update._prepare_database_staging(
                    active_root, clone_existing=True
                )
                try:
                    with update._using_database_path(stage):
                        update._build_database_artifacts(
                            {'demo': ['adk']}, {'demo': 'MLST'}
                        )
                    self.assertTrue(update.blast_database_is_ready(active_root, deep=True))
                    self.assertTrue(update.blast_database_is_ready(stage, deep=True))
                    self.assertEqual(
                        (active_root / 'mlst.fasta').read_bytes(), active_fasta
                    )
                finally:
                    shutil.rmtree(stage, ignore_errors=True)
            finally:
                update.pathdb = old_pathdb


if __name__ == '__main__':
    unittest.main()
