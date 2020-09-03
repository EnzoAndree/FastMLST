#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import logging
import xml.etree.ElementTree as ET
from Bio import SeqIO
from pickle import dump
from pickle import load
from datetime import date
from subprocess import PIPE
from sys import version_info
from subprocess import Popen
from collections import defaultdict
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool
from tqdm import tqdm  # pip3 install tqdm
try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve
if version_info[0] < 3:
    from pathlib2 import Path  # pip2 install pathlib2
else:
    from pathlib import Path

logger = logging.getLogger('update_mlst')

pathdb = Path(str(Path(__file__).parent.absolute()) + '/' + 'pubmlst')
necessary_file = ['mlst.fasta.nhr', 'mlst.fasta.nsq', 'mlst.fasta.nin',
                  'scheme_number.pkl']


def save_obj(obj, name):
    with open(name, 'wb') as f:
        dump(obj, f, 2)


def load_obj(name):
    with open(name, 'rb') as f:
        return load(f)


def parseXML(xml):
    tree = ET.parse(xml)
    root = tree.getroot()
    species = defaultdict()
    for parent in root.iter('species'):
        data = []
        for child in parent.iter('url'):
            data.append(child.text)
        codename = data[1].strip('/').split('/')[-4].split('_')[1]
        species[codename] = data[1:]
    return species


def download_fasta(files):
    codename = files[0].strip('/').split('/')[-4].split('_')[1]
    outdir = Path(str(pathdb) + '/schemes' + '/' + codename)
    outdir.mkdir(exist_ok=True, parents=True)
    for file in files:
        if 'profiles_csv' in file:
            urlretrieve(file, str(outdir) + '/' +
                        file.strip('/').split('/')[-4].split('_')[1] + '.txt')
        else:
            urlretrieve(file, str(outdir) + '/' +
                        file.strip('/').split('/')[-2] + '.tfa')
    return '[OK] {}'.format(codename)


def update_mlstdb(threads):
    pathdb.mkdir(exist_ok=True, parents=True)
    urlretrieve('https://pubmlst.org/data/dbases.xml', str(pathdb) +
                '/dbases.xml')
    logger.info('https://pubmlst.org/data/dbases.xml downloaded to {}'
                .format(str(pathdb) + '/dbases.xml'))
    datadb = parseXML(str(pathdb) + '/dbases.xml')
    logger.info('Starting download of all schemes')
    t = ThreadPool(threads)
    genome_mlst = []
    for result in tqdm(t.imap_unordered(download_fasta, datadb.values()),
                       total=len(datadb.values()),
                       desc='Downloading Schemes using {} threads'.
                       format(threads), unit='Schemes', leave=True):
        genome_mlst.append(result)
    t.close()
    t.join()
    logger.info('Schemes were downloaded')
    fastas = Path(str(pathdb) + '/schemes').glob('*/*.tfa')
    allfasta = []
    scheme_number = defaultdict()
    for species, data in datadb.items():
        scheme_number[species] = [genes.strip('/').split('/')[-2].
                                  split('.')[0] for genes in data[1:]]
    save_obj(scheme_number, str(pathdb) + '/scheme_number.pkl')
    logger.info('Schemes object was created in {}'.format(str(pathdb) +
                                                          '/scheme_number.pkl')
                )
    for fasta in fastas:
        for record in SeqIO.parse(str(fasta), 'fasta'):
            scheme = str(fasta.parent).split('/')[-1]
            record.id = '{}.{}'.format(scheme, record.id)
            record.description = ''
            allfasta.append(record)
    outfna = 'mlst.fasta'
    SeqIO.write(allfasta, str(pathdb) + '/' + outfna, 'fasta')
    blastdb_cmd = 'makeblastdb -hash_index -in {0} -dbtype nucl -title \
                   "PubMLST_{1}" -parse_seqids'
    blastdb_cmd = blastdb_cmd.format(str(pathdb) + '/' +
                                     outfna, date.today().strftime("%d%m%y"))
    DB_process = Popen(blastdb_cmd, shell=True, stdin=PIPE, stdout=PIPE,
                       stderr=PIPE)
    DB_process.wait()
    logger.info('BLASTdb was created using pubmlst data')
    logger.info('Update PubMLST Complete')
