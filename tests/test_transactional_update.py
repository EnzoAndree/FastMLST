import os
import shutil
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from fastmlst import update_mlst_kit as update


class FakeClient:
    def close(self):
        return None


class TransactionalDownloadTests(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.old_pathdb = update.pathdb
        self.old_credentials = update._oauth_credentials
        update.pathdb = Path(self.tmpdir) / 'pubmlst'

    def tearDown(self):
        update.pathdb = self.old_pathdb
        update._oauth_credentials = self.old_credentials
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    @staticmethod
    def item():
        return {
            'codename': 'pubmlst_demo_seqdef_1',
            'database': 'pubmlst_demo_seqdef',
            'scheme_id': 1,
            'profiles_csv': 'https://example.test/profiles',
            'loci': ['https://example.test/loci/adk'],
            'description': 'MLST',
            'remote_version': {
                'records': 1,
                'locus_count': 1,
                'last_updated': '2026-01-01',
                'last_added': '2026-01-01',
            },
        }

    def test_failed_scheme_download_preserves_installed_copy(self):
        item = self.item()
        scheme_dir = update.pathdb / 'schemes' / item['codename']
        scheme_dir.mkdir(parents=True)
        old_fasta = scheme_dir / 'adk.tfa'
        old_fasta.write_text('>old\nAAAA\n', encoding='utf-8')

        with patch.object(update, 'fetch_text', side_effect=RuntimeError('offline')):
            ok, *_rest = update.download_scheme_data((item, None))

        self.assertFalse(ok)
        self.assertEqual(old_fasta.read_text(encoding='utf-8'), '>old\nAAAA\n')
        self.assertFalse(list((update.pathdb / 'schemes').glob('.*.download-*')))

    def test_successful_scheme_download_removes_stale_loci(self):
        item = self.item()
        scheme_dir = update.pathdb / 'schemes' / item['codename']
        scheme_dir.mkdir(parents=True)
        (scheme_dir / 'obsolete.tfa').write_text('>1\nAAAA\n', encoding='utf-8')

        def fake_text(url, content_type='text/plain', session=None):
            if url.endswith('/profiles'):
                return 'ST\tadk\n1\t1\n'
            return '>1\nACGT\n'

        with patch.object(update, 'fetch_text', side_effect=fake_text):
            with patch.object(
                update,
                'fetch_json',
                return_value={
                    'id': 'adk',
                    'alleles_fasta': 'https://example.test/alleles/adk',
                },
            ):
                ok, *_rest = update.download_scheme_data((item, None))

        self.assertTrue(ok)
        self.assertTrue((scheme_dir / 'adk.tfa').is_file())
        self.assertFalse((scheme_dir / 'obsolete.tfa').exists())

    def test_incomplete_local_scheme_does_not_match_remote(self):
        item = self.item()
        item['loci'].append('https://example.test/loci/glyA')
        item['remote_version']['locus_count'] = 2
        scheme_dir = update.pathdb / 'schemes' / item['codename']
        scheme_dir.mkdir(parents=True)
        (scheme_dir / f"{item['codename']}.txt").write_text(
            'ST\tadk\tglyA\n1\t1\t1\n', encoding='utf-8'
        )
        (scheme_dir / 'adk.tfa').write_text('>1\nACGT\n', encoding='utf-8')
        update.save_local_scheme_meta(
            item['codename'],
            item['database'],
            item['scheme_id'],
            item['description'],
            item['remote_version'],
        )
        with patch.object(
            update, 'fetch_scheme_remote_version', return_value=item['remote_version']
        ):
            matches = update.local_scheme_matches_remote(
                item['codename'], item['database'], item['scheme_id']
            )
        self.assertFalse(matches)

    def test_scheme_validation_rejects_profile_missing_locus_column(self):
        item = self.item()
        scheme_dir = update.pathdb / 'staged-scheme'
        scheme_dir.mkdir(parents=True)
        (scheme_dir / f"{item['codename']}.txt").write_text(
            'ST\tglyA\n1\t1\n', encoding='utf-8'
        )
        (scheme_dir / 'adk.tfa').write_text('>adk_1\nACGT\n', encoding='utf-8')
        with self.assertRaisesRegex(RuntimeError, 'missing locus column'):
            update.validate_scheme_directory(scheme_dir, item, ['adk'])

    def test_transaction_failure_leaves_active_database_untouched(self):
        update.pathdb.mkdir(parents=True)
        marker = update.pathdb / 'dbases.xml'
        marker.write_text('old database', encoding='utf-8')
        update._oauth_credentials = {
            'client_id': 'id',
            'client_secret': 'secret',
            'auth_db': update.AUTH_TEST_DB,
        }

        with patch.object(update, '_build_authenticated_session', return_value=FakeClient()):
            with patch.object(update, 'build_scheme_catalog_direct', return_value=[self.item()]):
                with patch.object(
                    update,
                    '_update_database_in_place',
                    side_effect=RuntimeError('download failed'),
                ):
                    with self.assertRaisesRegex(RuntimeError, 'download failed'):
                        update.update_mlstdb_selected(1, selectors=['db:1'])

        self.assertEqual(marker.read_text(encoding='utf-8'), 'old database')
        parent = update.pathdb.parent
        self.assertFalse(list(parent.glob(f'.{update.pathdb.name}.staging-*')))
        self.assertFalse(list(parent.glob(f'.{update.pathdb.name}.backup-*')))

    def test_successful_transaction_replaces_complete_tree(self):
        update.pathdb.mkdir(parents=True)
        (update.pathdb / 'old.txt').write_text('old', encoding='utf-8')
        update._oauth_credentials = {
            'client_id': 'id',
            'client_secret': 'secret',
            'auth_db': update.AUTH_TEST_DB,
        }

        def fake_update(*_args, **_kwargs):
            (update.pathdb / 'new.txt').write_text('new', encoding='utf-8')

        with patch.object(update, '_build_authenticated_session', return_value=FakeClient()):
            with patch.object(update, 'build_scheme_catalog_direct', return_value=[self.item()]):
                with patch.object(update, '_update_database_in_place', side_effect=fake_update):
                    update.update_mlstdb_selected(1, selectors=['db:1'])

        self.assertEqual((update.pathdb / 'new.txt').read_text(encoding='utf-8'), 'new')
        self.assertTrue((update.pathdb / 'old.txt').is_file())

    def _seed_installed_scheme(self):
        scheme_dir = update.pathdb / 'schemes' / 'installed'
        scheme_dir.mkdir(parents=True)
        (scheme_dir / 'installed.txt').write_text(
            'ST\told\n1\t1\n', encoding='utf-8'
        )
        (scheme_dir / 'old.tfa').write_text('>1\nAACCGGTT\n', encoding='utf-8')
        update.save_obj({'installed': ['old']}, update.pathdb / 'scheme_number.pkl')
        update.save_obj({'installed': 'MLST'}, update.pathdb / update.SCHEME_LIST_FILE)
        (update.pathdb / 'sentinel.txt').write_text('active-old', encoding='utf-8')

    @staticmethod
    def _fake_download_into_current_stage(task):
        item, _session = task
        scheme_dir = update.pathdb / 'schemes' / item['codename']
        scheme_dir.mkdir(parents=True, exist_ok=True)
        (scheme_dir / f"{item['codename']}.txt").write_text(
            'ST\tadk\n1\t1\n', encoding='utf-8'
        )
        (scheme_dir / 'adk.tfa').write_text('>1\nACGTACGT\n', encoding='utf-8')
        return True, item['codename'], ['adk'], item['description'], None, None

    @staticmethod
    def _successful_makeblastdb(command, **_kwargs):
        if command[0] == 'blastdbcmd':
            return SimpleNamespace(returncode=0, stdout='Database: test', stderr='')
        prefix = Path(command[command.index('-out') + 1])
        for suffix in ('.nhr', '.nin', '.nsq'):
            Path(str(prefix) + suffix).write_text('index', encoding='utf-8')
        return SimpleNamespace(returncode=0, stdout='', stderr='')

    def test_real_transaction_pipeline_builds_and_publishes_consistent_artifacts(self):
        self._seed_installed_scheme()
        update._oauth_credentials = {
            'client_id': 'id',
            'client_secret': 'secret',
            'auth_db': update.AUTH_TEST_DB,
        }
        with patch.object(update, '_build_authenticated_session', return_value=FakeClient()):
            with patch.object(update, 'build_scheme_catalog_direct', return_value=[self.item()]):
                with patch.object(update, 'local_scheme_matches_remote', return_value=False):
                    with patch.object(
                        update, 'download_scheme_data', side_effect=self._fake_download_into_current_stage
                    ):
                        with patch.object(
                            update.subprocess,
                            'run',
                            side_effect=self._successful_makeblastdb,
                        ):
                            update.update_mlstdb_selected(1, selectors=['db:1'])

        manifest = update.load_obj(update.pathdb / 'scheme_number.pkl')
        self.assertEqual(manifest['installed'], ['old'])
        self.assertEqual(manifest[self.item()['codename']], ['adk'])
        self.assertTrue(update.blast_database_is_ready())
        self.assertEqual(
            update._count_fasta_records(update.pathdb / 'mlst.fasta'), 2
        )

    def test_makeblastdb_failure_rolls_back_entire_update(self):
        self._seed_installed_scheme()
        update._oauth_credentials = {
            'client_id': 'id',
            'client_secret': 'secret',
            'auth_db': update.AUTH_TEST_DB,
        }
        failed = SimpleNamespace(returncode=1, stdout='', stderr='bad database')
        with patch.object(update, '_build_authenticated_session', return_value=FakeClient()):
            with patch.object(update, 'build_scheme_catalog_direct', return_value=[self.item()]):
                with patch.object(update, 'local_scheme_matches_remote', return_value=False):
                    with patch.object(
                        update, 'download_scheme_data', side_effect=self._fake_download_into_current_stage
                    ):
                        with patch.object(update.subprocess, 'run', return_value=failed):
                            with self.assertRaisesRegex(RuntimeError, 'bad database'):
                                update.update_mlstdb_selected(1, selectors=['db:1'])

        self.assertEqual(
            (update.pathdb / 'sentinel.txt').read_text(encoding='utf-8'),
            'active-old',
        )
        self.assertFalse(
            (update.pathdb / 'schemes' / self.item()['codename']).exists()
        )

    def test_directory_swap_restores_old_tree_if_new_rename_fails(self):
        update.pathdb.mkdir(parents=True)
        (update.pathdb / 'old.txt').write_text('old', encoding='utf-8')
        stage = Path(tempfile.mkdtemp(
            prefix=f'.{update.pathdb.name}.staging-', dir=str(update.pathdb.parent)
        ))
        (stage / 'new.txt').write_text('new', encoding='utf-8')
        real_replace = update.os.replace

        def fail_new_publish(source, destination):
            if (
                Path(source).name == stage.name
                and Path(destination).name == update.pathdb.name
            ):
                raise OSError('rename failed')
            return real_replace(source, destination)

        with patch.object(update.os, 'replace', side_effect=fail_new_publish):
            with self.assertRaisesRegex(OSError, 'rename failed'):
                update._publish_database_staging(update.pathdb, stage)

        self.assertEqual(
            (update.pathdb / 'old.txt').read_text(encoding='utf-8'), 'old'
        )
        self.assertTrue(stage.is_dir())
        self.assertFalse(list(update.pathdb.parent.glob(f'.{update.pathdb.name}.backup-*')))

    def test_recovery_restores_backup_and_removes_stale_stage(self):
        backup = update.pathdb.parent / f'.{update.pathdb.name}.backup-test'
        backup.mkdir(parents=True)
        (backup / 'old.txt').write_text('old', encoding='utf-8')
        stale_stage = update.pathdb.parent / f'.{update.pathdb.name}.staging-test'
        stale_stage.mkdir()

        update._recover_interrupted_database_swap(update.pathdb)

        self.assertEqual(
            (update.pathdb / 'old.txt').read_text(encoding='utf-8'), 'old'
        )
        self.assertFalse(backup.exists())
        self.assertFalse(stale_stage.exists())

    def test_recovery_does_not_lose_backup_if_empty_target_was_recreated(self):
        update.pathdb.mkdir(parents=True)
        backup = update.pathdb.parent / f'.{update.pathdb.name}.backup-test'
        backup.mkdir(parents=True)
        (backup / 'old.txt').write_text('old', encoding='utf-8')

        update._recover_interrupted_database_swap(update.pathdb)

        self.assertEqual(
            (update.pathdb / 'old.txt').read_text(encoding='utf-8'), 'old'
        )
        self.assertFalse(backup.exists())

    def test_writer_cannot_swap_while_reader_holds_shared_lock(self):
        handle = update.acquire_database_read_lock(update.pathdb)
        try:
            with self.assertRaisesRegex(RuntimeError, 'already using'):
                with update._database_update_lock(update.pathdb):
                    self.fail('writer lock should not be acquired')
        finally:
            update.release_database_read_lock(handle)

        with update._database_update_lock(update.pathdb):
            pass

    def test_reader_lock_supports_read_only_shared_database(self):
        update.pathdb.mkdir(parents=True)
        parent = update.pathdb.parent
        lock_path = parent / f'.{update.pathdb.name}{update.UPDATE_LOCK_SUFFIX}'
        lock_path.write_text('shared database lock\n', encoding='utf-8')
        parent_mode = parent.stat().st_mode
        lock_mode = lock_path.stat().st_mode
        handle = None
        try:
            lock_path.chmod(0o444)
            parent.chmod(0o555)
            handle = update.acquire_database_read_lock(update.pathdb)
            self.assertIsNotNone(handle)
        finally:
            update.release_database_read_lock(handle)
            parent.chmod(parent_mode)
            lock_path.chmod(lock_mode)

    def test_reader_fails_safely_when_shared_lock_cannot_be_created(self):
        with patch.object(
            Path,
            'open',
            side_effect=[FileNotFoundError(), PermissionError('read-only parent')],
        ):
            with self.assertRaisesRegex(RuntimeError, 'administrator to create'):
                update.acquire_database_read_lock(update.pathdb)

    def test_writer_lock_is_readable_despite_restrictive_umask(self):
        previous_umask = os.umask(0o077)
        try:
            with update._database_update_lock(update.pathdb):
                pass
        finally:
            os.umask(previous_umask)
        lock_path = (
            update.pathdb.parent
            / f'.{update.pathdb.name}{update.UPDATE_LOCK_SUFFIX}'
        )
        self.assertEqual(lock_path.stat().st_mode & 0o444, 0o444)


class CatalogPolicyTests(unittest.TestCase):
    def test_catalog_auth_failure_names_database_registration_to_fix(self):
        root = [{
            'databases': [{
                'name': 'pubmlst_new_seqdef',
                'href': 'https://rest.pubmlst.org/db/pubmlst_new_seqdef',
            }]
        }]
        with patch.object(
            update,
            'fetch_json',
            side_effect=[
                root,
                update.PubMLSTAuthenticationError('HTTP 401'),
            ],
        ):
            with self.assertRaises(RuntimeError) as caught:
                update.fetch_public_scheme_catalog(session='authenticated')

        message = str(caught.exception)
        self.assertIn('pubmlst_new_seqdef', message)
        self.assertIn('may mean', message)
        self.assertIn('MY ACCOUNT > Database registrations', message)
        self.assertIn('fastmlst --pubmlst-reset-auth', message)
        self.assertIn('existing catalog was not replaced', message)

    def test_catalog_discovery_rejects_missing_schemes_payload(self):
        root = [{
            'databases': [{
                'name': 'pubmlst_demo_seqdef',
                'href': 'https://rest.pubmlst.org/db/pubmlst_demo_seqdef',
            }]
        }]
        with patch.object(update, 'fetch_json', side_effect=[root, {}]):
            with self.assertRaisesRegex(RuntimeError, 'catalog discovery was incomplete'):
                update.fetch_public_scheme_catalog(session='authenticated')

    def test_catalog_discovery_rejects_scheme_without_url(self):
        root = [{
            'databases': [{
                'name': 'pubmlst_demo_seqdef',
                'href': 'https://rest.pubmlst.org/db/pubmlst_demo_seqdef',
            }]
        }]
        schemes = {'records': 1, 'schemes': [{'description': 'MLST'}]}
        with patch.object(update, 'fetch_json', side_effect=[root, schemes]):
            with self.assertRaisesRegex(RuntimeError, 'scheme entry has no scheme URL'):
                update.fetch_public_scheme_catalog(session='authenticated')

    def test_missing_catalog_is_fetched_without_writing_active_tree(self):
        catalog = [{'codename': 'fetched'}]
        with patch.object(update, 'load_scheme_catalog', return_value=None):
            with patch.object(
                update, 'fetch_public_scheme_catalog', return_value=catalog
            ) as fetch:
                with patch.object(update, 'save_scheme_catalog') as save:
                    result, from_cache = update._load_or_build_scheme_catalog(
                        session='authenticated'
                    )
        self.assertEqual(result, catalog)
        self.assertFalse(from_cache)
        fetch.assert_called_once_with(session='authenticated')
        save.assert_not_called()

    def test_full_update_excludes_large_schemes_unless_requested(self):
        standard = {
            'description': 'MLST',
            'remote_version': {'locus_count': 7},
            'loci': list(range(7)),
        }
        large = {
            'description': 'cgMLST',
            'remote_version': {'locus_count': 2000},
            'loci': list(range(2000)),
        }
        selected, excluded = update._filter_full_catalog([standard, large])
        self.assertEqual(selected, [standard])
        self.assertEqual(excluded, [large])
        selected_all, excluded_all = update._filter_full_catalog(
            [standard, large], include_large_schemes=True
        )
        self.assertEqual(selected_all, [standard, large])
        self.assertEqual(excluded_all, [])


if __name__ == '__main__':
    unittest.main()
