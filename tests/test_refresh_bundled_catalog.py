import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from fastmlst import update_mlst_kit


SCRIPT_PATH = Path(__file__).resolve().parents[1] / 'scripts' / 'refresh_bundled_catalog.py'
SPEC = importlib.util.spec_from_file_location('refresh_bundled_catalog', SCRIPT_PATH)
refresh_bundled_catalog = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(refresh_bundled_catalog)


class TestRefreshBundledCatalog(unittest.TestCase):
    def test_requires_configured_authentication_before_building_session(self):
        with patch.object(
            update_mlst_kit,
            'configure_pubmlst_auth_from_environment',
            side_effect=RuntimeError('credentials are missing'),
        ) as configure_auth:
            with patch.object(update_mlst_kit, '_build_authenticated_session') as build_session:
                with patch.object(update_mlst_kit, 'fetch_public_scheme_catalog') as fetch_catalog:
                    with self.assertRaisesRegex(SystemExit, 'without PubMLST authentication'):
                        refresh_bundled_catalog.main()

        configure_auth.assert_called_once_with(require_auth=True)
        build_session.assert_not_called()
        fetch_catalog.assert_not_called()

    def test_rejects_missing_authenticated_session(self):
        with patch.object(
            update_mlst_kit,
            'configure_pubmlst_auth_from_environment',
        ) as configure_auth:
            with patch.object(
                update_mlst_kit,
                '_build_authenticated_session',
                return_value=None,
            ):
                with patch.object(update_mlst_kit, 'fetch_public_scheme_catalog') as fetch_catalog:
                    with self.assertRaisesRegex(SystemExit, 'authenticated PubMLST session'):
                        refresh_bundled_catalog.main()

        configure_auth.assert_called_once_with(require_auth=True)
        fetch_catalog.assert_not_called()

    def test_fetches_and_saves_catalog_with_authenticated_session(self):
        old_pathdb = update_mlst_kit.pathdb
        catalog = [{'codename': 'test_scheme'}]
        session = object()

        with tempfile.TemporaryDirectory() as tmpdir:
            fake_module = Path(tmpdir) / 'fastmlst' / 'update_mlst_kit.py'
            with patch.object(update_mlst_kit, '__file__', str(fake_module)):
                with patch.object(
                    update_mlst_kit,
                    'configure_pubmlst_auth_from_environment',
                ) as configure_auth:
                    with patch.object(
                        update_mlst_kit,
                        '_build_authenticated_session',
                        return_value=session,
                    ):
                        with patch.object(
                            update_mlst_kit,
                            'fetch_public_scheme_catalog',
                            return_value=catalog,
                        ) as fetch_catalog:
                            with patch.object(update_mlst_kit, 'save_scheme_catalog') as save_catalog:
                                refresh_bundled_catalog.main()

            expected_bundle = fake_module.parent / 'bundle'
            self.assertTrue(expected_bundle.is_dir())

        configure_auth.assert_called_once_with(require_auth=True)
        fetch_catalog.assert_called_once_with(session=session)
        save_catalog.assert_called_once_with(catalog)
        self.assertEqual(update_mlst_kit.pathdb, old_pathdb)


if __name__ == '__main__':
    unittest.main()
