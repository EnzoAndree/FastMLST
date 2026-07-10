#!/usr/bin/env python

import argparse
import logging
import multiprocessing
import os
from codecs import decode
from collections import defaultdict
from itertools import repeat
from multiprocessing import Pool
from multiprocessing import cpu_count
from pathlib import Path
from re import compile
from sys import exit
from sys import stderr
from sys import stdout

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm  # pip3 install tqdm

from fastmlst import __version__
from fastmlst.mlst import MLST
from fastmlst.update_mlst_kit import necessary_file
import fastmlst.update_mlst_kit as update_mlst_kit


def unescaped_str(arg_str):
    return decode(str(arg_str), 'unicode_escape')


def check_coverage_range(value):
    fvalue = float(value)
    if fvalue <= 0 or fvalue > 100:
        raise argparse.ArgumentTypeError(
            'The coverage velue must to be [0..100]'
        )
    return fvalue


def check_identity_range(value):
    fvalue = float(value)
    if fvalue <= 0 or fvalue > 100:
        raise argparse.ArgumentTypeError(
            'The identity velue must to be [0..100]'
        )
    return fvalue


def check_threads(value):
    ivalue = int(value)
    if ivalue < 1:
        raise argparse.ArgumentTypeError('The number of threads must be >= 1')
    return ivalue


def safe_reset_database(db_path, logger):
    from shutil import copy2, rmtree
    import tempfile

    db_path = Path(db_path).resolve()
    if not db_path.exists():
        return
    if not db_path.is_dir():
        raise RuntimeError(f'Database path is not a directory: {db_path}')
    if db_path == Path('/'):
        raise RuntimeError('Refusing to remove root directory "/"')

    markers = {
        'dbases.xml',
        'schemes',
        'mlst.fasta',
        'mlst.fasta.nhr',
        'mlst.fasta.nsq',
        'mlst.fasta.nin',
        update_mlst_kit.SCHEME_CATALOG_FILE,
        update_mlst_kit.SCHEME_CATALOG_META_FILE,
    }
    current_entries = {p.name for p in db_path.iterdir()}
    is_empty = len(current_entries) == 0
    has_markers = bool(current_entries.intersection(markers))
    if not is_empty and not has_markers:
        raise RuntimeError(
            f'Refusing to remove non-empty path without FastMLST markers: {db_path}'
        )

    catalog_names = (
        update_mlst_kit.SCHEME_CATALOG_FILE,
        update_mlst_kit.SCHEME_CATALOG_META_FILE,
    )
    tmp = tempfile.mkdtemp()
    preserved = []
    try:
        for fname in catalog_names:
            src = db_path / fname
            if src.is_file():
                copy2(src, Path(tmp) / fname)
                preserved.append(fname)

        logger.info(f'Removing existing FastMLST database directory: {db_path}')
        rmtree(str(db_path))

        if preserved:
            db_path.mkdir(parents=True, exist_ok=True)
            for fname in preserved:
                copy2(Path(tmp) / fname, db_path / fname)
            logger.info(
                'Kept PubMLST scheme catalog cache (%s) for faster --update-mlst',
                ', '.join(preserved),
            )
    finally:
        rmtree(tmp, ignore_errors=True)


def runMLST(margument):
    genome, cov, ident, sep, header, shcheme = margument
    return MLST(genome, cov, ident, sep, header, shcheme)


def resolve_update_selectors_arg(update_mlst_arg, legacy_select_schemes_arg=''):
    if update_mlst_arg not in (None, ''):
        return update_mlst_arg
    return legacy_select_schemes_arg


def main():
    V = f'%(prog)s v{__version__}'
    parser = argparse.ArgumentParser(
        description='⚡️🧬 FastMLST: A multi-core tool for multilocus sequence typing of draft genome assemblies'
    )
    parser.add_argument(type=str, nargs='*', dest='genomes')
    parser.add_argument('-t', '--threads', type=check_threads, default=cpu_count(),
                        help='Number of threads to use (default {})'.
                        format(cpu_count()))
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        choices=[0, 1, 2],
                        help='Verbose output level choices: [0, 1, 2]')
    parser.add_argument('-s', '--separator', type=unescaped_str, default=',',
                        help='Choose a character to use as a separator' +
                        ' (default ",")')
    parser.add_argument('-sch', '--scheme', type=str,
                        help='Set a scheme target (I am not dumb, let me choose a scheme by myself!)')
    parser.add_argument('--scheme-list', action='store_true',
                        help='List PubMLST schemes (cache, bundled catalog, or one-time API fetch if missing)')
    parser.add_argument('--scheme-list-update', action='store_true',
                        help='Refresh scheme catalog from the PubMLST API, then print the list')
    parser.add_argument('--list-remote-schemes', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--list-remote-schemes-update', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--installed-scheme-stats', action='store_true',
                        help='List each installed scheme on one line: PubMLST description, type, ST count, alleles, per-locus counts')
    parser.add_argument('-fo', '--fastaoutput', type=str, default='',
                        help='File name of the concatenated alleles output' +
                        ' (default "")')
    parser.add_argument('-to', '--tableoutput', type=str, default=stdout,
                        help='File name of the MLST table output' +
                        ' (default STDOUT)')
    parser.add_argument('-cov', '--coverage', type=check_coverage_range,
                        default=99,
                        help='DNA %%Cov to report high quality partial allele [?]' +
                        ' (default 99%%)')
    parser.add_argument('-pid', '--identity', type=check_identity_range,
                        default=95,
                        help='DNA %%Identity of full allelle to consider' +
                        ' \'similar\' [~] (default 95%%)')
    parser.add_argument('--update-mlst', type=str, nargs='?', const='',
                        help='Update PubMLST from the API. Pass ALL for the full catalog (very slow), '
                             'or a comma-separated list of database:scheme_id / stable codenames.')
    parser.add_argument('--select-schemes', type=str, default='',
                        help=argparse.SUPPRESS)
    parser.add_argument('--redownload-all-schemes', action='store_true',
                        help='Ignore local version metadata and re-download every scheme (default: skip schemes that match remote API metadata)')
    parser.add_argument('-sp', '--splited-output', type=str, default='',
                        help='Directory output for splited alleles' +
                        ' (default "")')
    parser.add_argument('--fasta2line', action='store_true',
                        help='The fasta files will be in fasta2line format')
    parser.add_argument('--longheader', action='store_true',
                        help='If --longheader is invoked, the header of FASTA' +
                        ' file contain a long description')
    parser.add_argument('--legacy', action='store_true',
                        help='If --legacy is invoked, the csv reported contain the gene name' +
                        ' and the allele id in the row [adk(1),atpA(4),dxr(7),glyA(1),recA(1),sodA(3),tpi(3)].' +
                        ' This option is only available when the --scheme is defined')
    parser.add_argument('-n', '--novel', type=str,
                        help='File name of the novel alleles')
    parser.add_argument('-V', '--version', action='version',
                        version=V, help='Show program\'s version number and exit')
    parser.add_argument(
        '--db_path',
        type=str,
        default=None,
        help='Custom directory for MLST database (default: ~/.cache/fastmlst/pubmlst)'
    )
    parser.add_argument('--pubmlst-client-id', type=str, default=None,
                        help='PubMLST OAuth client ID (or FASTMLST_PUBMLST_CLIENT_ID)')
    parser.add_argument('--pubmlst-client-secret', type=str, default=None,
                        help='PubMLST OAuth client secret (or FASTMLST_PUBMLST_CLIENT_SECRET)')
    parser.add_argument('--pubmlst-connect', action='store_true',
                        help='Run one-time PubMLST OAuth setup and save credentials/tokens')
    args = parser.parse_args()

    if args.db_path:
        update_mlst_kit.set_pathdb(args.db_path)

    client_id = args.pubmlst_client_id or os.environ.get('FASTMLST_PUBMLST_CLIENT_ID')
    client_secret = args.pubmlst_client_secret or os.environ.get('FASTMLST_PUBMLST_CLIENT_SECRET')
    access_token = os.environ.get('FASTMLST_PUBMLST_ACCESS_TOKEN')
    access_secret = os.environ.get('FASTMLST_PUBMLST_ACCESS_SECRET')
    verifier = os.environ.get('FASTMLST_PUBMLST_VERIFIER')
    update_mlst_kit.set_oauth_credentials(client_id, client_secret, access_token, access_secret, verifier)
    if client_id and client_secret:
        update_mlst_kit.save_pubmlst_client_credentials(client_id, client_secret)
    if args.pubmlst_connect:
        update_mlst_kit.connect_pubmlst(client_id=client_id, client_secret=client_secret, verifier=verifier)
        logger = logging.getLogger('FastMLST')
        logger.info('PubMLST OAuth credentials/tokens are configured.')
        exit()

    selector_arg = resolve_update_selectors_arg(args.update_mlst, args.select_schemes)
    selectors = [s.strip() for s in selector_arg.split(',') if s.strip()]
    scheme_list_refresh = (
        args.scheme_list_update or args.list_remote_schemes_update
    )
    if (
        args.scheme_list
        or scheme_list_refresh
        or args.list_remote_schemes
    ):
        has_cached_catalog = update_mlst_kit.load_scheme_catalog() is not None
        session = None
        if scheme_list_refresh or not has_cached_catalog:
            session = update_mlst_kit._build_authenticated_session()
        catalog, from_cache, used_live_api = update_mlst_kit.get_remote_scheme_catalog(
            session=session,
            force_refresh=scheme_list_refresh,
        )
        if scheme_list_refresh:
            if not used_live_api:
                print(
                    'WARNING: --scheme-list-update did not contact the live PubMLST API '
                    '(using existing catalog data only).',
                    file=stderr,
                )
            elif session is None:
                print(
                    'WARNING: No PubMLST OAuth session; refresh uses anonymous API access. '
                    'Some databases may return HTTP 401 — run --pubmlst-connect for full access.',
                    file=stderr,
                )
        meta = update_mlst_kit.load_scheme_catalog_meta()
        update_mlst_kit.print_remote_scheme_catalog(
            catalog,
            from_cache=from_cache,
            meta=meta,
        )
        exit()

    if args.installed_scheme_stats:
        update_mlst_kit.print_installed_scheme_stats()
        exit()

    split_namefromcode = compile(r'(?P<gene>.+)\((?P<novel>~?)(?P<number>\d+)(?P<partial>\??)\)')

    if args.verbose == 0:
        logging.basicConfig(level=logging.WARNING,
                            format='[%(asctime)s] %(levelname)s@%(name)s: %(message)s',
                            datefmt='%H:%M:%S')
        logger = logging.getLogger('FastMLST')
    elif args.verbose == 1:
        logging.basicConfig(level=logging.INFO,
                            format='[%(asctime)s] %(levelname)s@%(name)s: %(message)s',
                            datefmt='%H:%M:%S')
        logger = logging.getLogger('FastMLST')
    elif args.verbose == 2:
        logging.basicConfig(level=logging.DEBUG,
                            format='[%(asctime)s] %(levelname)s@%(name)s: %(message)s',
                            datefmt='%H:%M:%S')
        logger = logging.getLogger('FastMLST')
    update_mlst_kit.pathdb.mkdir(exist_ok=True, parents=True)
    is_all_files = all((update_mlst_kit.pathdb / f).is_file() for f in necessary_file)
    if not is_all_files and args.update_mlst is None:
        print(
            'ERROR: PubMLST database files are missing under {}.\n'
            'Run: fastmlst --update-mlst ALL\n'
            '  (full catalog, very slow) or e.g.\n'
            '  fastmlst --update-mlst "pubmlst_cdifficile_seqdef:1"'.format(
                update_mlst_kit.pathdb
            ),
            file=stderr,
        )
        exit(2)

    if args.update_mlst is not None:
        from fastmlst.update_mlst_kit import parse_update_mlst_selectors
        from fastmlst.update_mlst_kit import update_mlstdb_selected
        try:
            mlst_download_selectors = parse_update_mlst_selectors(selectors)
        except ValueError as exc:
            parser.error(str(exc))
        if mlst_download_selectors is None:
            safe_reset_database(update_mlst_kit.pathdb, logger)
        update_mlstdb_selected(
            args.threads,
            selectors=mlst_download_selectors,
            force_redownload_schemes=args.redownload_all_schemes,
        )
        exit()
    if not args.genomes:
        parser.print_help(stderr)
        exit()
    if args.scheme is not None:
        args.scheme = args.scheme.lower()
        scheme_dir = update_mlst_kit.pathdb / 'schemes'
        if args.scheme in [d.name for d in scheme_dir.iterdir()]:
            logger.info('Ok my little buddy, i trust your judgment. I will ' +
                        f'proceed with the search using only the following scheme: {args.scheme}')
        else:
            logger.error(f'Are you sure that "{args.scheme}" is a supported scheme?')
            logger.error('Don\'t worry my little buddy. You are a human ' +
                         'after all. I\'ll keep trying to choose the best scheme.')
            args.scheme = None
    genome_mlst = []
    multipleargs = list(zip(args.genomes,
                            repeat(args.coverage),
                            repeat(args.identity),
                            repeat(args.separator),
                            repeat(args.longheader),
                            repeat(args.scheme),
                            ))
    with Pool(args.threads) as p:
        for result in tqdm(p.imap(runMLST, multipleargs),
                           total=len(multipleargs),
                           desc='Scanning Genomes using {} threads'.
                           format(args.threads), unit='Genomes', leave=False):
            genome_mlst.append(result)
    fastaconcat = []
    fastasplited = defaultdict(list)
    fastanovelconcat = []
    str_alleles = ''
    dict_alleles = []
    for genome in genome_mlst:
        if genome.blastresult:
            if not genome.descarted \
               and not genome.contamination \
               and not genome.allelemissing:
                fastaconcat.append(genome.concat_alleles)
                if args.splited_output != '':
                    for allele in genome.name_alleles:
                        fasta = genome.alleles[allele]
                        fasta.id = genome.beautiname
                        fasta.description = f'{allele}({genome.dict_st[allele]})'
                        fastasplited[allele].append(fasta)

                if args.novel and genome.novel_alleles:
                    for novelallele in genome.novel_alleles:
                        genenovel = split_namefromcode.search(novelallele)
                        if genenovel:
                            gene_name = genenovel.group('gene')
                            try:
                                seq = genome.alleles[gene_name]
                            except KeyError as e:
                                logger.warning(f"Novel allele '{gene_name}' not found in extracted alleles for {genome.beautiname}: {e}")
                                continue
                        else:
                            logger.warning(f"No gene match found for novel allele entry: {novelallele}")
                            continue
                        seq.id = novelallele + '@' + genome.beautiname
                        seq.description = ''
                        fastanovelconcat.append(seq)
            str_alleles += genome.str_st
            str_alleles += '\n'
            dict_alleles.append(genome.dict_st)
    if args.fasta2line and args.fastaoutput != '':
        SeqIO.write(fastaconcat, args.fastaoutput, 'fasta-2line')
        if args.novel:
            SeqIO.write(fastanovelconcat, args.novel, 'fasta-2line')
    elif args.fastaoutput != '':
        SeqIO.write(fastaconcat, args.fastaoutput, 'fasta')
        if args.novel:
            SeqIO.write(fastanovelconcat, args.novel, 'fasta')
    if args.fasta2line and args.splited_output != '':
        spout = Path(args.splited_output)
        spout.mkdir(exist_ok=True, parents=True)
        for gene, fastalist in fastasplited.items():
            SeqIO.write(fastalist, f'{spout.absolute()}/{gene}.fasta', 'fasta-2line')
    elif args.splited_output != '':
        spout = Path(args.splited_output)
        spout.mkdir(exist_ok=True, parents=True)
        for gene, fastalist in fastasplited.items():
            SeqIO.write(fastalist, f'{spout.absolute()}/{gene}.fasta', 'fasta')
    if type(args.tableoutput) == str:
        if args.scheme is not None:
            if args.legacy:
                with open(args.tableoutput, 'w') as output_handle:
                    print(str_alleles[:-1], file=output_handle)
            else:
                df = pd.DataFrame(dict_alleles)
                df.to_csv(f'{args.tableoutput}', index=False, sep=args.separator)
        else:
            with open(args.tableoutput, 'w') as output_handle:
                print(str_alleles[:-1], file=output_handle)
    else:
        if args.scheme is not None:
            if args.legacy:
                print(str_alleles[:-1], file=args.tableoutput)
            else:
                df = pd.DataFrame(dict_alleles)
                args.tableoutput.write(df.to_csv(index=False, sep=args.separator))
        else:
            print(str_alleles[:-1], file=args.tableoutput)


if __name__ == '__main__':
    multiprocessing.set_start_method('fork')
    main()
