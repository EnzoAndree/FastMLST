import shutil
import tempfile
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from fastmlst.cli import resolve_update_selectors_arg
from fastmlst.cli import safe_reset_database
from fastmlst import mlst as mlst_module
from fastmlst.mlst import MLST
from fastmlst import update_mlst_kit


class DummyLogger:
    def info(self, *_args, **_kwargs):
        return None


class TestMlstBugFixes(unittest.TestCase):
    def test_allele_number_matches_uses_exact_tokens(self):
        self.assertTrue(MLST.allele_number_matches("1", "1"))
        self.assertTrue(MLST.allele_number_matches("1", "~1"))
        self.assertTrue(MLST.allele_number_matches("1", "1?"))
        self.assertTrue(MLST.allele_number_matches("10", "1|10"))
        self.assertFalse(MLST.allele_number_matches("1", "10"))
        self.assertFalse(MLST.allele_number_matches("2", "12|20"))

    def test_stassignment_returns_ambiguous_sentinel(self):
        tmpdir = tempfile.mkdtemp()
        old_pathdb = update_mlst_kit.pathdb
        try:
            scheme_dir = Path(tmpdir) / "schemes" / "testscheme"
            scheme_dir.mkdir(parents=True, exist_ok=True)
            profile = scheme_dir / "testscheme.txt"
            profile.write_text("ST\tadk\n1\t1\n2\t1\n", encoding="utf-8")
            update_mlst_kit.set_pathdb(tmpdir)

            obj = MLST.__new__(MLST)
            obj.scheme = "testscheme"
            obj.score = {"scheme": {"adk": "1"}}
            obj.blastn_cli = "blastn"
            obj.fasta = "sample.fna"
            result = obj.STassignment()
            self.assertEqual(result, "ambiguous_ST")
        finally:
            update_mlst_kit.set_pathdb(old_pathdb)
            shutil.rmtree(tmpdir)

    def test_make_blast_sets_failure_state_on_non_zero_exit(self):
        obj = MLST.__new__(MLST)
        obj.fasta_opened = ">x\nACGT\n"
        obj.fasta = "sample.fna"
        obj.blastresult = True

        with patch("fastmlst.mlst.subprocess.run") as mocked_run:
            mocked_run.return_value.returncode = 1
            mocked_run.return_value.stdout = ""
            mocked_run.return_value.stderr = "blast failed"
            result = obj.make_blast()

        self.assertIsNone(result)
        self.assertFalse(obj.blastresult)


class TestCliSafetyFixes(unittest.TestCase):
    def test_update_mlst_uses_inline_selector_value(self):
        self.assertEqual(
            resolve_update_selectors_arg("ALL", "pubmlst_cdifficile_seqdef_1"),
            "ALL",
        )

    def test_update_mlst_falls_back_to_legacy_select_schemes(self):
        self.assertEqual(
            resolve_update_selectors_arg("", "pubmlst_cdifficile_seqdef_1"),
            "pubmlst_cdifficile_seqdef_1",
        )

    def test_safe_reset_database_rejects_non_marker_directory(self):
        tmpdir = tempfile.mkdtemp()
        try:
            p = Path(tmpdir) / "custom"
            p.mkdir(parents=True, exist_ok=True)
            (p / "notes.txt").write_text("do not delete", encoding="utf-8")
            with self.assertRaises(RuntimeError):
                safe_reset_database(p, DummyLogger())
            self.assertTrue(p.exists())
        finally:
            shutil.rmtree(tmpdir)

    def test_safe_reset_database_allows_marker_directory(self):
        tmpdir = tempfile.mkdtemp()
        try:
            p = Path(tmpdir) / "pubmlst"
            p.mkdir(parents=True, exist_ok=True)
            (p / "dbases.xml").write_text("<root/>", encoding="utf-8")
            safe_reset_database(p, DummyLogger())
            self.assertFalse(p.exists())
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_safe_reset_database_preserves_scheme_catalog(self):
        tmpdir = tempfile.mkdtemp()
        try:
            p = Path(tmpdir) / "pubmlst"
            p.mkdir(parents=True, exist_ok=True)
            (p / "dbases.xml").write_text("<root/>", encoding="utf-8")
            (p / update_mlst_kit.SCHEME_CATALOG_FILE).write_text(
                '[{"codename": "x"}]', encoding="utf-8"
            )
            (p / update_mlst_kit.SCHEME_CATALOG_META_FILE).write_text(
                '{"count": 1}', encoding="utf-8"
            )
            safe_reset_database(p, DummyLogger())
            self.assertTrue(p.is_dir())
            self.assertFalse((p / "dbases.xml").exists())
            self.assertTrue((p / update_mlst_kit.SCHEME_CATALOG_FILE).is_file())
            self.assertTrue((p / update_mlst_kit.SCHEME_CATALOG_META_FILE).is_file())
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_table_output_uses_given_handle_for_dataframe_mode(self):
        table_buf = StringIO()
        df = __import__('pandas').DataFrame([{'Genome': 'g1', 'ST': 1}])
        table_buf.write(df.to_csv(index=False, sep=','))
        self.assertEqual(table_buf.getvalue(), 'Genome,ST\ng1,1\n')


class TestDynamicPathdb(unittest.TestCase):
    def test_set_pathdb_updates_mlst_runtime_lookup(self):
        old_pathdb = update_mlst_kit.pathdb
        tmpdir = tempfile.mkdtemp()
        try:
            update_mlst_kit.set_pathdb(tmpdir)
            self.assertEqual(mlst_module.update_mlst_kit.pathdb, Path(tmpdir))
        finally:
            update_mlst_kit.set_pathdb(old_pathdb)
            shutil.rmtree(tmpdir, ignore_errors=True)


class TestDownloadSchemeDataHandling(unittest.TestCase):
    def test_download_scheme_data_empty_codename_returns_failure(self):
        item = {'codename': '  ', 'database': 'pubmlst_x_seqdef', 'scheme_id': 1}
        ok, *_rest = update_mlst_kit.download_scheme_data((item, None))
        self.assertFalse(ok)

    def test_download_scheme_data_no_alleles_fasta_does_not_raise(self):
        tmpdir = tempfile.mkdtemp()
        old_pathdb = update_mlst_kit.pathdb
        try:
            update_mlst_kit.pathdb = Path(tmpdir)
            item = {
                'codename': 'pubmlst_testfail_seqdef_9',
                'database': 'pubmlst_testfail_seqdef',
                'scheme_id': 9,
                'profiles_csv': 'http://example.invalid/profiles',
                'loci': ['https://rest.pubmlst.org/db/pubmlst_testfail_seqdef/loci/1'],
                'description': 'test',
            }
            with patch.object(update_mlst_kit, 'fetch_text', return_value='ST\tadk\n1\t1\n'):
                with patch.object(
                    update_mlst_kit,
                    'fetch_json',
                    return_value={'id': 'adk', 'alleles_fasta': None},
                ):
                    ok, *_r = update_mlst_kit.download_scheme_data((item, None))
            self.assertFalse(ok)
            scheme_dir = Path(tmpdir) / 'schemes' / item['codename']
            self.assertFalse(scheme_dir.is_dir())
        finally:
            update_mlst_kit.pathdb = old_pathdb
            shutil.rmtree(tmpdir, ignore_errors=True)


class TestInstalledSchemeStats(unittest.TestCase):
    def test_print_installed_scheme_stats_one_line(self):
        tmpdir = tempfile.mkdtemp()
        old_pathdb = update_mlst_kit.pathdb
        try:
            root = Path(tmpdir)
            scheme_dir = root / 'schemes' / 'pubmlst_foo_seqdef_1'
            scheme_dir.mkdir(parents=True)
            (scheme_dir / 'pubmlst_foo_seqdef_1.txt').write_text(
                'ST\tadk\tglyA\n1\t1\t2\n2\t2\t1\n',
                encoding='utf-8',
            )
            (scheme_dir / 'adk.tfa').write_text('>1\nAT\n>2\nGC\n', encoding='utf-8')
            (scheme_dir / 'glyA.tfa').write_text('>1\nAA\n', encoding='utf-8')
            update_mlst_kit.pathdb = root
            update_mlst_kit.save_obj(
                {'pubmlst_foo_seqdef_1': 'Clostridioides difficile MLST (7 locus)'},
                str(root / update_mlst_kit.SCHEME_LIST_FILE),
            )
            from io import StringIO

            buf = StringIO()
            with patch('sys.stdout', buf):
                update_mlst_kit.print_installed_scheme_stats()
            out = buf.getvalue().strip()
            self.assertIn('pubmlst_foo_seqdef_1', out)
            self.assertIn('type=MLST', out)
            self.assertIn(
                'description=Clostridioides difficile MLST (7 locus)',
                out,
            )
            self.assertIn('STs=2', out)
            self.assertIn('alleles=3', out)
            self.assertIn('adk=2', out)
            self.assertIn('glyA=1', out)
            self.assertEqual(out.count('\n'), 0)
        finally:
            update_mlst_kit.pathdb = old_pathdb
            shutil.rmtree(tmpdir, ignore_errors=True)


class TestRmlstExcludedFromConcat(unittest.TestCase):
    def test_scheme_excluded_from_mlst_concat(self):
        self.assertTrue(
            update_mlst_kit.scheme_excluded_from_mlst_concat('pubmlst_rmlst_seqdef_1')
        )
        self.assertFalse(
            update_mlst_kit.scheme_excluded_from_mlst_concat('pubmlst_cdifficile_seqdef_1')
        )


class TestParseUpdateMlstSelectors(unittest.TestCase):
    def test_all_returns_none(self):
        self.assertIsNone(update_mlst_kit.parse_update_mlst_selectors(['ALL']))
        self.assertIsNone(update_mlst_kit.parse_update_mlst_selectors(['all']))

    def test_specific_schemes_passthrough(self):
        s = ['pubmlst_cdifficile_seqdef:1']
        self.assertEqual(update_mlst_kit.parse_update_mlst_selectors(s), s)

    def test_empty_raises(self):
        with self.assertRaises(ValueError):
            update_mlst_kit.parse_update_mlst_selectors([])

    def test_all_mixed_raises(self):
        with self.assertRaises(ValueError):
            update_mlst_kit.parse_update_mlst_selectors(['ALL', 'pubmlst_x_seqdef:1'])


class TestLoadOrBuildSchemeCatalog(unittest.TestCase):
    def test_load_or_build_uses_cache_without_fetch(self):
        sentinel = [{'codename': 'cached'}]
        with patch.object(update_mlst_kit, 'load_scheme_catalog', return_value=sentinel):
            with patch.object(
                update_mlst_kit,
                'fetch_public_scheme_catalog',
            ) as fetch_mock:
                cat, used = update_mlst_kit._load_or_build_scheme_catalog(session=None)
        fetch_mock.assert_not_called()
        self.assertIs(cat, sentinel)
        self.assertTrue(used)


class TestGetRemoteSchemeCatalog(unittest.TestCase):
    def test_force_refresh_always_calls_fetch(self):
        sentinel = [{'codename': 'x', 'database': 'd', 'scheme_id': 1}]
        with patch.object(
            update_mlst_kit,
            'fetch_public_scheme_catalog',
            return_value=sentinel,
        ) as fetch_mock:
            with patch.object(update_mlst_kit, 'save_scheme_catalog') as save_mock:
                cat, from_cache, used_api = update_mlst_kit.get_remote_scheme_catalog(
                    session=None,
                    force_refresh=True,
                )
        fetch_mock.assert_called_once()
        save_mock.assert_called_once()
        self.assertIs(cat, sentinel)
        self.assertFalse(from_cache)
        self.assertTrue(used_api)

    def test_without_refresh_uses_cache_and_skips_fetch(self):
        sentinel = [{'codename': 'cached'}]
        with patch.object(
            update_mlst_kit,
            'load_scheme_catalog',
            return_value=sentinel,
        ):
            with patch.object(
                update_mlst_kit,
                'fetch_public_scheme_catalog',
            ) as fetch_mock:
                cat, from_cache, used_api = update_mlst_kit.get_remote_scheme_catalog(
                    force_refresh=False,
                )
        fetch_mock.assert_not_called()
        self.assertIs(cat, sentinel)
        self.assertTrue(from_cache)
        self.assertFalse(used_api)


class TestBundledSchemeCatalog(unittest.TestCase):
    def test_load_scheme_catalog_falls_back_to_package_bundle(self):
        tmpdir = tempfile.mkdtemp()
        old_pathdb = update_mlst_kit.pathdb
        try:
            update_mlst_kit.pathdb = Path(tmpdir)
            bundled = (
                update_mlst_kit._package_bundle_dir()
                / update_mlst_kit.SCHEME_CATALOG_FILE
            )
            if not bundled.is_file():
                self.skipTest('package bundle/scheme_catalog.json not present')
            cat = update_mlst_kit.load_scheme_catalog()
            self.assertIsInstance(cat, list)
            self.assertGreater(len(cat), 0)
            self.assertIn('codename', cat[0])
            self.assertFalse(
                (Path(tmpdir) / update_mlst_kit.SCHEME_CATALOG_FILE).is_file()
            )
        finally:
            update_mlst_kit.pathdb = old_pathdb
            shutil.rmtree(tmpdir, ignore_errors=True)


class TestDirectSchemeSelector(unittest.TestCase):
    def test_parse_direct_scheme_selector(self):
        db, sid = update_mlst_kit.parse_direct_scheme_selector('pubmlst_cdifficile_seqdef:1')
        self.assertEqual(db, 'pubmlst_cdifficile_seqdef')
        self.assertEqual(sid, 1)
        db2, sid2 = update_mlst_kit.parse_direct_scheme_selector(
            'pubmlst_cdifficile_seqdef_1'
        )
        self.assertEqual(db2, 'pubmlst_cdifficile_seqdef')
        self.assertEqual(sid2, 1)
        self.assertEqual(update_mlst_kit.parse_direct_scheme_selector('bad')[0], None)
        self.assertEqual(update_mlst_kit.parse_direct_scheme_selector('foo:bar')[0], None)
        self.assertEqual(
            update_mlst_kit.parse_direct_scheme_selector('pubmlst_cdifficile_seqdef')[0],
            None,
        )

    def test_scheme_type_from_description(self):
        self.assertEqual(
            update_mlst_kit.scheme_type_from_description('Clostridioides difficile MLST'),
            'MLST',
        )
        self.assertEqual(
            update_mlst_kit.scheme_type_from_description('Salmonella cgMLST scheme'),
            'cgMLST',
        )
        self.assertEqual(update_mlst_kit.scheme_type_from_description(''), 'unknown')

    def test_ensure_scheme_codename_never_empty(self):
        code = update_mlst_kit.ensure_scheme_codename(
            'pubmlst_cdifficile_seqdef', 1, ''
        )
        self.assertEqual(code, 'pubmlst_cdifficile_seqdef_1')

    def test_ensure_scheme_codename_ignores_generic_description(self):
        """Descriptions like 'Extended MLST' must not become a shared folder name."""
        code = update_mlst_kit.ensure_scheme_codename(
            'pubmlst_cbotulinum_seqdef', 2, 'Extended MLST'
        )
        self.assertEqual(code, 'pubmlst_cbotulinum_seqdef_2')

    def test_remote_version_from_scheme_detail(self):
        v = update_mlst_kit.remote_version_from_scheme_detail({
            'id': 1,
            'records': 1308,
            'locus_count': 7,
            'last_updated': '2025-03-01',
        })
        self.assertEqual(v['records'], 1308)
        self.assertEqual(v['locus_count'], 7)
        self.assertEqual(v['last_updated'], '2025-03-01')


class TestOAuthConfigFixes(unittest.TestCase):
    def test_set_oauth_credentials_requires_all_fields(self):
        tmpdir = tempfile.mkdtemp()
        old_dir = update_mlst_kit.CONFIG_DIR
        old_file = update_mlst_kit.CONFIG_FILE
        try:
            update_mlst_kit.CONFIG_DIR = Path(tmpdir)
            update_mlst_kit.CONFIG_FILE = Path(tmpdir) / "pubmlst_oauth.json"
            update_mlst_kit.set_oauth_credentials()
            with self.assertRaises(ValueError):
                update_mlst_kit.set_oauth_credentials(client_id="abc")
            with self.assertRaises(ValueError):
                update_mlst_kit.set_oauth_credentials(client_id="a", client_secret="b", access_token="tok")
        finally:
            update_mlst_kit.CONFIG_DIR = old_dir
            update_mlst_kit.CONFIG_FILE = old_file
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_save_and_load_pubmlst_credentials(self):
        tmpdir = tempfile.mkdtemp()
        old_dir = update_mlst_kit.CONFIG_DIR
        old_file = update_mlst_kit.CONFIG_FILE
        try:
            update_mlst_kit.CONFIG_DIR = Path(tmpdir)
            update_mlst_kit.CONFIG_FILE = Path(tmpdir) / "pubmlst_oauth.json"
            update_mlst_kit.save_pubmlst_client_credentials("cid", "csecret")
            cid, csecret = update_mlst_kit.load_pubmlst_client_credentials()
            self.assertEqual(cid, "cid")
            self.assertEqual(csecret, "csecret")
        finally:
            update_mlst_kit.CONFIG_DIR = old_dir
            update_mlst_kit.CONFIG_FILE = old_file
            shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
