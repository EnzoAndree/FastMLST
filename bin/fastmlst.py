#!/usr/bin/env python

import argparse
import logging
from sys import exit
from sys import stderr
from sys import stdout
from Bio import SeqIO
from multiprocessing import Pool
from multiprocessing import cpu_count
from fastmlst.update_mlst_kit import pathdb
from fastmlst.update_mlst_kit import load_obj
from fastmlst.update_mlst_kit import necessary_file
from itertools import repeat
from fastmlst.mlst import MLST
from codecs import decode
from tqdm import tqdm  # pip3 install tqdm
from pathlib import Path
import pandas as pd
from re import compile
from collections import defaultdict
import multiprocessing
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


def runMLST(margument):
    genome, cov, ident, sep, header, shcheme = margument
    return MLST(genome, cov, ident, sep, header, shcheme)


def main():
    V = '%(prog)s v0.0.19'
    parser = argparse.ArgumentParser(
        description='⚡️🧬 FastMLST: A multi-core tool for multilocus sequence typing of draft genome assemblies'
    )
    parser.add_argument(type=str, nargs='*', dest='genomes')
    parser.add_argument('-t', '--threads', type=int, default=cpu_count(),
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
                        help='Show all schemes supported')
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
    parser.add_argument('--update-mlst', action='store_true',
                        help='Perform an update of the PubMLST database')
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
    args = parser.parse_args()

    # If the user provided a custom database path, update it
    if args.db_path:
        update_mlst_kit.set_pathdb(args.db_path)

    # Verbose?
    formatter = logging.Formatter('[%(asctime)s] %(levelname)s@%(name)s: %(message)s')
    ch = logging.StreamHandler()
    split_namefromcode = compile(r'(?P<gene>.+)\((?P<novel>~?)(?P<number>\d+)(?P<partial>\??)\)') # OMG genename can be alphanumeric

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
    # Check for pubmlst
    update_mlst_kit.pathdb.mkdir(exist_ok=True, parents=True)
    is_all_files = all((update_mlst_kit.pathdb / f).is_file() for f in necessary_file)
    # If update_mlst is true o any necesary file are missing update pubmlst
    if args.update_mlst or not is_all_files:
        from shutil import rmtree
        from fastmlst.update_mlst_kit import update_mlstdb
        rmtree(str(update_mlst_kit.pathdb))
        update_mlstdb(args.threads)
        if args.update_mlst:
            exit()
    if args.scheme_list:
        from fastmlst.update_mlst_kit import show_scheme_list
        show_scheme_list()
        exit()
    if not args.genomes:
        parser.print_help(stderr)
        exit()
    # Check if there are a target scheme
    if args.scheme != None:
        args.scheme = args.scheme.lower()
        scheme_dir = update_mlst_kit.pathdb/'schemes'
        if args.scheme in [d.name for d in scheme_dir.iterdir()]:
            logger.info('Ok my little buddy, i trust your judgment. I will '+
                          f'proceed with the search using only the following scheme: {args.scheme}')
        else:
            logger.error(f'Are you sure that "{args.scheme}" is a supported scheme?')
            logger.error('Don\'t worry my little buddy. You are a human '+
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
    # Show the shcheme list and exit
    fastaconcat = []
    fastasplited = defaultdict(list)
    fastanovelconcat = []
    str_alleles = ''
    dict_alleles = []
    for genome in genome_mlst:
        # only export if
        # 1. Allele not contain Ns
        # 2. No alleles missing
        # 3. No contamination in genome
        # if genome.blastresult\
        #         and not genome.descarted\
        #         and not genome.contamination\
        #         and not genome.allelemissing:
        if genome.blastresult:
            if not genome.descarted\
               and not genome.contamination\
               and not genome.allelemissing:
                fastaconcat.append(genome.concat_alleles)
                # splited alleles
                if args.splited_output != '':
                    for allele in genome.name_alleles:
                        fasta = genome.alleles[allele]
                        fasta.id = genome.beautiname
                        fasta.description = f'{allele}({genome.dict_st[allele]})'
                        if allele not in fastasplited.keys():
                            fastasplited[allele].append(fasta)
                        else:
                            fastasplited[allele].append(fasta)

                if args.novel and genome.novel_alleles:
                    for novelallele in genome.novel_alleles:
                        genenovel = split_namefromcode.search(novelallele)
                        if genenovel:
                            gene_name = genenovel.group('gene')
                            try:
                                # Debugging output
                                seq = genome.alleles[gene_name]
                            except KeyError as e:
                                print(f"KeyError: {e} - The key '{gene_name}' was not found in genome.alleles.")
                                continue  # Skip this allele and continue with the next
                        else:
                            print(f"No match found for novelallele: {novelallele}")
                            continue
                        seq.id = novelallele + '@' + genome.beautiname
                        seq.description = ''
                        fastanovelconcat.append(seq)
            str_alleles += genome.str_st
            str_alleles += '\n'
            dict_alleles.append(genome.dict_st)
    # FastMLSTv0.0.12 by default do not write the fasta file 
    if args.fasta2line and args.fastaoutput != '':
        SeqIO.write(fastaconcat, args.fastaoutput, 'fasta-2line')
        if args.novel:
            SeqIO.write(fastanovelconcat, args.novel, 'fasta-2line')
    elif args.fastaoutput != '':
        SeqIO.write(fastaconcat, args.fastaoutput, 'fasta')
        if args.novel:
            SeqIO.write(fastanovelconcat, args.novel, 'fasta')
    # Novel output in version v0.0.14
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
    # output formated
    if type(args.tableoutput) == str:
        if args.scheme != None:
            if args.legacy:
                print(str_alleles[:-1], file=open(args.tableoutput, 'w'))
            else:
                df = pd.DataFrame(dict_alleles)
                df.to_csv(f'{args.tableoutput}', index=False, sep=args.separator)
        else:
            print(str_alleles[:-1], file=open(args.tableoutput, 'w'))
    else:
        if args.scheme != None:
            if args.legacy:
                print(str_alleles[:-1], file=args.tableoutput)
            else:
                df = pd.DataFrame(dict_alleles)
                print(df.to_csv(index=False, sep=args.separator))
        else:
            print(str_alleles[:-1], file=args.tableoutput)

if __name__ == '__main__':
    multiprocessing.set_start_method('fork')  # or 'spawn' or 'forkserver'
    main()
