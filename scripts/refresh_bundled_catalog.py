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
    session = u._build_authenticated_session()
    if session:
        print('Using saved PubMLST OAuth session for catalog fetch.', flush=True)
    else:
        print('Using anonymous API access (may fail with HTTP 401 on some databases).', flush=True)
    print('Fetching catalog from https://rest.pubmlst.org (this may take several minutes)...', flush=True)
    catalog = u.fetch_public_scheme_catalog(session=session)
    prev = u.pathdb
    try:
        u.pathdb = bundle
        u.save_scheme_catalog(catalog)
    finally:
        u.pathdb = prev
    print(f'Wrote {len(catalog)} schemes to {bundle}', flush=True)


if __name__ == '__main__':
    main()
