import logging
import xml.etree.ElementTree as ET
from Bio import SeqIO
from pickle import dump
from pickle import load
from datetime import date
from subprocess import PIPE
from subprocess import Popen
from collections import defaultdict
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool
from tqdm import tqdm  # pip3 install tqdm
from urllib.request import urlretrieve
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

def best_guess_codename(text):
    # i hope that this is correct
    text = text.replace('.', '').strip().lower()
    text = text.replace('candidatus ', '')
    text = text.replace('/', '') # the dots are used in blast_filter to separate the scheme form data, use other sep
    text = text.replace('(', '')
    text = text.replace(')', '')
    # no more bugs plz!
    if len(text.split(' ')) == 1:
        # is a genus scheme
        if '#' in text:
            # sh#t, there is multiples schemes for a sigle specie :(
            components = text.split('#')
            genus = components[0]
            schemenumber = components[-1]
            codename = f'{genus}#{schemenumber}'
        else:
            codename = f'{text}'
    elif 'spp' in text:
        # is a genus scheme
        if '#' in text:
            # sh#t, there is multiples schemes for a sigle specie :(
            components = text.split(' ')
            genus = components[0]
            schemenumber = text.split('#')[-1]
            codename = f'{genus}#{schemenumber}'
        else:
            components = text.split(' ')
            genus = components[0]
            codename = f'{genus}'
    elif len(text.split(' ')) == 2:
        # is a clasic genus_species like scheme
        if '#' in text:
            # sh#t, there is multiples schemes for a sigle specie :(
            components = text.split(' ')
            genus = components[0][0]
            specie = components[1].split('#')[0]
            schemenumber = text.split('#')[-1]
            codename = f'{genus}{specie}#{schemenumber}'
        else:
            components = text.split(' ')
            genus = components[0][0]
            specie = components[1]
            codename = f'{genus}{specie}'
    elif len(text.split(' ')) > 2:
        # is a genus_species_etc... like scheme
        if '#' in text:
            # sh#t, there is multiples schemes for a sigle specie :(
            components = text.split(' ')
            genus = components[0][0]
            specie = components[1].split('#')[0]
            extra = '_'.join(components[2:])
            schemenumber = text.split('#')[-1]
            codename = f'{genus}{specie}_{extra}#{schemenumber}'
        else:
            components = text.split(' ')
            genus = components[0][0]
            specie = components[1]
            extra = '_'.join(components[2:])
            codename = f'{genus}{specie}_{extra}'
    return codename.strip()


def parseXML(xml):
    tree = ET.parse(xml)
    root = tree.getroot()
    species = defaultdict()
    for parent in root.iter('species'):
        data = []
        for child in parent.iter('url'):
            data.append(child.text)
        # BUG solved!, this is not the scheme code calculate it using best_guess_codename
        # codename = '_'.join(data[1].strip('/').split('/')[-4].split('_')[1:-1])
        codename = best_guess_codename(parent.text)
        species[codename] = data[1:]
    return species


def download_fasta(items):
    codename = items[0]
    files = items[1]
    outdir = Path(str(pathdb) + '/schemes' + '/' + codename)
    outdir.mkdir(exist_ok=True, parents=True)
    for file in files:
        if 'profiles_csv' in file:
            out_filename = codename
            urlretrieve(file, str(outdir) + '/' +
                        out_filename + '.txt')
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
    with ThreadPool(threads) as t:
        for result in tqdm(t.imap(download_fasta, datadb.items()),
                           total=len(datadb.items()),
                           desc='Downloading Schemes using {} threads'.
                           format(threads), unit='Schemes', leave=True):
            genome_mlst.append(result)
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
                                     outfna, date.today().strftime('%d%m%y'))
    DB_process = Popen(blastdb_cmd, shell=True, stdin=PIPE, stdout=PIPE,
                       stderr=PIPE)
    DB_process.wait()
    logger.info('BLASTdb was created using pubmlst data')
    logger.info('Update PubMLST Complete')

def show_scheme_list():
    datadb = Path(str(pathdb) + '/dbases.xml')
    if not datadb.is_file():
        # update the database
        from sys import exit
        logger.error('There is no dbases.xml, please update the database')
        exit()
    tree = ET.parse(datadb)
    root = tree.getroot()
    species = defaultdict()
    for parent in root.iter('species'):
        data = []
        for child in parent.iter('url'):
            data.append(child.text)
        # BUG solved!, this is not the scheme code calculate it using best_guess_codename
        # codename = '_'.join(data[1].strip('/').split('/')[-4].split('_')[1:-1])
        codename = best_guess_codename(parent.text)
        species[codename] = parent.text
    print(f'There are {len(species)} schemes (A round of applause to @keithajolley! (Jolley, et al., 2018)):\n')
    i = 1
    for sch, species in species.items():
        print(f'({i}) {sch}: {species.strip()}')
        i += 1
