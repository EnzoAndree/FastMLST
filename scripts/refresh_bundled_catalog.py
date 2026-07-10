#!/usr/bin/env python3
"""
Regenerate fastmlst/bundle/scheme_catalog.json and scheme_catalog_meta.json
from the live PubMLST REST API. Run from the repository root (may take several minutes).

  python scripts/refresh_bundled_catalog.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import fastmlst.update_mlst_kit as u  # noqa: E402


def main() -> None:
    bundle = Path(u.__file__).resolve().parent / 'bundle'
    bundle.mkdir(parents=True, exist_ok=True)

    try:
        u.configure_pubmlst_auth_from_environment(require_auth=True)
    except (RuntimeError, ValueError) as exc:
        raise SystemExit(
            f'Cannot refresh the bundled catalog without PubMLST authentication: {exc}'
        ) from exc

    try:
        session = u._build_authenticated_session()
    except (RuntimeError, ValueError, OSError) as exc:
        raise SystemExit(f'Could not create an authenticated PubMLST session: {exc}') from exc
    if session is None:
        raise SystemExit(
            'Cannot refresh the bundled catalog without an authenticated PubMLST session.'
        )

    print('Using PubMLST OAuth session for catalog fetch.', flush=True)
    print('Fetching catalog from https://rest.pubmlst.org (this may take several minutes)...', flush=True)
    try:
        catalog = u.fetch_public_scheme_catalog(session=session)
        prev = u.pathdb
        try:
            u.pathdb = bundle
            u.save_scheme_catalog(catalog)
        finally:
            u.pathdb = prev
    finally:
        close = getattr(session, 'close', None)
        if callable(close):
            close()
    print(f'Wrote {len(catalog)} schemes to {bundle}', flush=True)


if __name__ == '__main__':
    main()
