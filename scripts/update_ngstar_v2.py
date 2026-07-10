#!/usr/bin/env python3

import argparse
import logging
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fastmlst.ngstar import default_ngstar_database_path
from fastmlst.ngstar import update_ngstar_v2_database


def main():
    parser = argparse.ArgumentParser(
        description='Download and build the independent NG-STAR v2 database.'
    )
    parser.add_argument(
        '--db-path',
        default=str(default_ngstar_database_path()),
        help='NG-STAR v2 database directory.',
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    update_ngstar_v2_database(args.db_path)


if __name__ == '__main__':
    main()
