from operator import itemgetter
import pandas as pd
from fastmlst.update_mlst_kit import pathdb
import logging
from fastmlst.update_mlst_kit import load_obj
from collections import defaultdict
from sys import exit
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from gzip import open as gopen
from io import StringIO
from pathlib import Path
import gzip
import bz2
import gc
import subprocess

magic_dict = {
    b"\x1f\x8b\x08": (gzip.open, 'rb'),
    b"\x42\x5a\x68": (bz2.BZ2File, 'r'),
    }

max_len = max(len(x) for x in magic_dict)

def open_by_magic(filename):
    with open(filename,  "rb") as f:
        file_start = f.read(max_len)
    for magic, (fn, flag) in magic_dict.items():
        if file_start.startswith(magic):
            return fn(filename, flag)
    return open(filename, 'r')


logger = logging.getLogger('mlst')

# Will be arbitrarily excluded. by Enzo.
excluded_by_default = ['abaumannii#2', 'ecoli#2', 'leptospira#2', 'leptospira#3',
                       'mgallisepticum#1', 'pmultocida#1', 'sthermophilus',
                       'vcholerae#2']

class MLST(object):
    def __init__(self, fasta, coverage=75, identity=95, sep=',',
            longheader=False, target_scheme=None):
        super(MLST, self).__init__()
        self.target_scheme = target_scheme
        self.longheader = longheader
        self.fasta = fasta
        self.fasta_opened = open_by_magic(self.fasta).read()
        if type(self.fasta_opened) != str:
            self.fasta_opened = self.fasta_opened.decode()
        self.coverage = coverage / 100.0
        self.identity = identity / 100.0
        self.sep = sep
        self.scheme_number = load_obj(str(pathdb) + '/scheme_number.pkl')
        self.beautiname = self.fasta.strip('/').split('/')[-1]
        self.blastn_cli = None
        # QCflags
        self.descarted = False
        self.contamination = False
        self.allelemissing = False
        self.blastresult = False
        # QCflags
        self.blast = self.make_blast()
        if self.blastresult:
            self.scheme = None
            self.score = None
            self.novel_alleles = []
            self.scoring()
            self.QCflags()
            if not self.allelemissing and\
                    not self.novel_alleles and\
                    not self.contamination:
                self.ST = self.STassignment()
            # if novel alleles, is new ST by default
            elif not self.allelemissing and\
                    self.novel_alleles and\
                    not self.contamination:
                self.ST = 'new_alleles'
            else:
                self.ST = '-'
            self.name_alleles = self.scheme_number[self.scheme]
            self.number_alleles = len(self.name_alleles)
            self.STnumber = None
            self.alleles = None
            self.concat_alleles = self.mlstex()
            if self.descarted:
                # If any allele has Ns or is broken in 2 contigs, do not determine STs
                self.ST = '-'
            self.str_st, self.dict_st = self.str_allelic_profile()

        # Release Ram!
        del self.fasta_opened
        del self.scheme_number
        del self.blast
        gc.collect()

    def __repr__(self,):
        return '{}–ST: {}'.format(self.beautiname, self.STnumber)

    def QCflags(self, ):
        for locus, value in self.score['scheme'].items():
            if '-' in value:
                self.allelemissing = True
            # mlst.py must check for contamination in mlstex()
            # elif '|' in value:
            #     # Check if is a duplication of same reported alleles
            #     if '~' in value:
            #         self.contamination = True
            #     else:
            #         valuelist = value.split('|')
            #         if len(set(valuelist)) == 1:
            #             pass
            #         else:
            #             self.contamination = True
            elif '~' in value or '?' in value:
                self.novel_alleles.append(f'{locus}({value})')


    def make_blast(self,):
        # Build the BLAST command without using the deprecated Bio.Blast.Application wrappers
        cmd = [
            "blastn",
            "-db", str(pathdb) + '/mlst.fasta',
            "-dust", "no",
            "-outfmt", "6 sseqid slen sstrand sstart send length nident gaps qseqid qstart qend",
            "-max_target_seqs", "130000",
            "-evalue", "1E-20",
            "-ungapped"
        ]
        self.blastn_cli = " ".join(cmd)
        logger.debug(self.blastn_cli + ' < ' + self.fasta)
        result = subprocess.run(cmd, input=self.fasta_opened, text=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = result.stdout, result.stderr
        if out == '':
            logger.warning(f'There is no result for: {self.blastn_cli} < {self.fasta}')
            return None
        self.blastresult = True
        blastfiltred = self.blast_filter(out)
        del out
        del err
        return blastfiltred

    def blast_filter(self, blast_out):
        header = ['sseqid', 'slen', 'sstrand', 'sstart', 'send', 'length',
                  'nident', 'gaps', 'qseqid', 'qstart', 'qend']
        dfblast = pd.read_csv(StringIO(blast_out), sep='\t', names=header)
        toint = ['slen', 'sstart', 'send', 'length', 'nident', 'gaps',
                 'qstart', 'qend']
        dfblast['coverage'] = (dfblast.length - dfblast.gaps) / dfblast.slen
        dfblast['identity'] = (dfblast.nident - dfblast.gaps) / dfblast.slen # this is a 'global' %identity
        dfblast[toint] = dfblast[toint].astype(int)
        # Bug Fixed on FastMLSTv0.0.12: The contigs identifier of some assemblers like unicycler is just a number, I opted to standardize all searches using strings
        dfblast['qseqid'] = dfblast['qseqid'].astype(str)
        # don't even look at the hits below these cov's and identities
        dfblast = dfblast.loc[(dfblast['coverage'] <= 1) & (dfblast['coverage'] >= self.coverage) & (dfblast['identity'] >= self.identity)] # insertions can not be processed properly yet
        # dfblast = dfblast.loc[dfblast['slen'] >= dfblast['length']] # if have an insertion slen < length
        if len(dfblast) == 0:
            # there is no result
            logger.warning(f'There is no result for: {self.blastn_cli} < {self.fasta}')
            self.blastresult = False
            return None
        else:
            dfblast = dfblast.join(
                dfblast['sseqid'].str.split('.', n=1, expand=True).
                rename(columns={0: 'scheme', 1: 'genenumber'}))
            dfblast = dfblast.join(
                dfblast['genenumber'].str.rsplit('_', n=1, expand=True).
                rename(columns={0: 'gene', 1: 'number'}))
            dfblast = dfblast.drop(['sseqid', 'genenumber'], axis=1)
            dfblast['genome_id'] = self.beautiname
            # dfblast.index = dfblast['genome_id']
            dfblast = dfblast[['genome_id', 'scheme', 'gene', 'number', 'slen',
                               'sstrand', 'sstart', 'send', 'length', 'nident',
                               'gaps', 'coverage', 'identity', 'qseqid',
                               'qstart', 'qend']]
            if self.target_scheme != None:
                dfblast = dfblast.loc[dfblast['scheme'] == self.target_scheme]
            else:
                dfblast = dfblast[~dfblast['scheme'].isin(excluded_by_default)] # decapitation!
            if len(dfblast) == 0:
                # there is no result
                logger.warning(f'There is no result for: {self.blastn_cli} < {self.fasta}')
                self.blastresult = False
                return None

            # dfblast.sort_index(inplace=True)
            # Grup by gene and select the best hit (cov=100% high ID)
            # Better Timing
            # dfblast = dfblast.sort_values(by=['coverage', 'nident', 'gaps' ],
            #       ascending=[False, False, True]).drop_duplicates(['gene'], keep='first')

            # genegrup = dfblast.groupby('gene')
            # blastfiltred_bygene = []
            # for gene, df_group in genegrup:
            #     df_group.sort_values(by=['coverage', 'nident', 'gaps' ],
            #                          ascending=[False, False, True], inplace=True)
            #     blastfiltred_bygene.append(df_group.head(5))
            # dfblast = pd.concat(blastfiltred_bygene, ignore_index=True)
            # del blastfiltred_bygene
            # del genegrup
            return dfblast


    def str_allelic_profile(self, ):
        if not isinstance(self.ST, pd.DataFrame):
            output = '{0}{3}{1}{3}{2}{3}'.format(self.beautiname,
                                                 self.scheme, self.ST,
                                                 self.sep)
            self.STnumber = self.ST
            dictofrows = {'Genome': self.beautiname,
                          'Scheme': self.scheme,
                          'ST': self.ST}
            for i in self.score['scheme'].keys():
                out = '{0}({1}){2}'.format(i, self.score['scheme'][i],
                                           self.sep)
                dictofrows[i] = self.score['scheme'][i]
                output += out
            output = output.strip(self.sep)
        else:
            output = '{0}{3}{1}{3}{2}{3}'.format(self.beautiname,
                                                 self.scheme,
                                                 self.ST.index.values[0],
                                                 self.sep)
            self.STnumber = self.ST.index.values[0]
            dictofrows = {'Genome': self.beautiname,
                          'Scheme': self.scheme,
                          'ST': self.ST.index.values[0]}
            for i in self.ST:
                out = '{0}({1}){2}'.format(i, self.ST[i].values[0],
                                           self.sep)
                dictofrows[i] = self.ST[i].values[0]
                output += out
            output = output.strip(self.sep)
        return (output, dictofrows)

    def is_context_complete(self, length, start, end):
        if start < 0 or end < 0:
            return False
        elif start > length or end > length:
            return False
        else:
            return True

    def STassignment(self, ):
        scheme_dir = str(pathdb) + '/schemes' + '/' + self.scheme
        STlist = Path(str(scheme_dir) + '/' + self.scheme + '.txt')
        dfSTlist = pd.read_csv(str(STlist), sep='\t', index_col=0)
        for key, value in self.score['scheme'].items():
            if '|' in value:
                value = list(set(value.split('|')))
                if len(value) == 1:
                    dfSTlist = dfSTlist.loc[dfSTlist[key] == int(value[0])]
                else:
                    return '-'
            else:
                dfSTlist = dfSTlist.loc[dfSTlist[key] == int(value)]
        if len(dfSTlist) == 1:
            return dfSTlist
        elif len(dfSTlist) == 0:
            return 'new_ST'
        else:
            logger.error('If you got here, congratulations, ' +
                         ' you found a place in maintenance STassignment()!')
            logger.error(dfSTlist, self.blastn_cli + ' < ' + self.fasta)
            logger.error(self.score['scheme'])

    def mlstex(self, ):
        fasta_output = dict()
        for record in SeqIO.parse(StringIO(self.fasta_opened), 'fasta'):
            try:
                pd_blast = self.blast.loc[(self.blast['qseqid'] == record.id) &
                                          (self.blast['scheme'] == self.scheme)
                                          ]
            except KeyError:
                continue
            if isinstance(pd_blast, pd.DataFrame):
                for row in pd_blast.iterrows():
                    if row[1]['number'] not in\
                            self.score['scheme'][row[1]['gene']] or\
                            row[1]['coverage'] < self.coverage:
                        continue
                    if row[1]['sstrand'] == 'plus':
                        if row[1]['slen'] == row[1]['send']:
                            # finish well
                            finishmissing = 0
                        else:
                            # is missing some nucleotides
                            finishmissing = row[1]['slen'] - row[1]['send']
                        if row[1]['sstart'] == 1:
                            # start well
                            startmissing = 0
                        else:
                            # is missing some nucleotides
                            startmissing = row[1]['sstart'] - 1
                        start = int(row[1]['qstart']) - 1 - startmissing
                        end = int(row[1]['qend']) + finishmissing
                        seq = record.seq[start:end]
                    else:
                        if row[1]['slen'] == row[1]['sstart']:
                            # start well
                            startmissing = 0
                        else:
                            # is missing some nucleotides
                            startmissing = row[1]['slen'] - row[1]['sstart']
                        if row[1]['send'] == 1:
                            # finish well
                            finishmissing = 0
                        else:
                            # is missing some nucleotides
                            finishmissing = row[1]['send'] - 1
                        start = int(row[1]['qstart']) - 1 - startmissing
                        end = int(row[1]['qend']) + finishmissing
                        seq = record.seq[start:end].reverse_complement()
                    if seq.count('N') > 0 or not \
                            self.is_context_complete(len(record), start, end):
                        self.descarted = True
                        continue
                    identificator = '{}|{}|{}|{}_{}'.format(self.beautiname,
                                                            record.id,
                                                            row[1]['gene'],
                                                            start, end)
                    record_fasta = SeqRecord(seq, id=identificator)
                    # BUG: SNP at end of aligment is detected as deletion. Probably 'fixed' in 'Grup by gene and select the best hit' commentary
                    # Check for equal gene name and diferent sequence (e.g. conaminations)
                    if row[1]['gene'] in fasta_output:
                        # gene name already in the output
                        if fasta_output[row[1]['gene']].seq == record_fasta.seq:
                            # Ok, it is the same allele
                            pass
                        else:
                            self.contamination = True
                    else:
                        fasta_output[row[1]['gene']] = record_fasta
                del pd_blast
            elif isinstance(pd_blast, pd.Series):
                logger.error('If you got here, congratulations, ' +
                             ' you found a place in maintenance mlstex()!')
        self.alleles = fasta_output
        if self.longheader:
            header = self.beautiname + ' '
            concatenatedseq = ''
            for genename in self.name_alleles:
                if genename in self.alleles.keys():
                    header += genename + '_'
                    concatenatedseq += self.alleles[genename].seq
            record_out = SeqRecord(concatenatedseq, id=header.strip('_'),
                                   description='Concatenated Sequences of MLST ' +
                                   'from ' + self.beautiname)
        else:
            concatenatedseq = ''
            for genename in self.name_alleles:
                if genename in self.alleles.keys():
                    concatenatedseq += self.alleles[genename].seq
            record_out = SeqRecord(concatenatedseq, id=self.beautiname,
                description='')
        return record_out

    def scoring(self, ):
        genome_query = set(self.blast['genome_id'].tolist())
        if len(genome_query) == 1:
            genome_query = list(genome_query)[0]
        else:
            logger.warning('Warning, more than one genome as query ',
                           genome_query)
            exit()
        # check completeness, perfect identity, snp identity,
        rank_list = defaultdict(dict)
        for scheme, group in self.blast.groupby('scheme'):
            rank_list[scheme] = defaultdict(dict)
            rank_list[scheme]['score'] = 0
            rank_list[scheme]['scheme'] = defaultdict()
            loci = self.scheme_number[scheme]
            N = len(loci)
            blast_scheme = self.blast[(self.blast.gene.isin(loci)) &
                                      (self.blast.scheme == scheme)].copy()
            blast_scheme.sort_values(by=['coverage', 'identity', 'length'],
                                     ascending=[False, False, False], inplace=True)
            for locus in loci:
                row = blast_scheme[blast_scheme.gene == locus]
                if len(row) == 0:
                    # allele missing
                    rank_list[scheme]['score'] += 0
                    rank_list[scheme]['scheme'][locus] = '-'
                elif len(row) == 1:
                    # only one allele
                    # if have an insertion slen < length
                    # if row['slen'].values[0] < row['length'].values[0]:
                    #     rank_list[scheme]['score'] += 20.0 / N
                    #     rank_list[scheme]['scheme'][locus] = \
                    #         '{}?'.format(row['number'].values[0])
                    if row['coverage'].values[0] == 1 and\
                            row['identity'].values[0] == 1:
                        # perfect match
                        rank_list[scheme]['score'] += 100.0 / N
                        rank_list[scheme]['scheme'][locus] = \
                            row['number'].values[0]
                    # BUG: Blast no make a full aligmnent, if a snp aries at the very end, the coverage is not 1
                    # In the fasta output, this is fixed lookat the input sequence, but this no update the blast table
                    elif row['coverage'].values[0] == 1 and\
                            row['identity'].values[0] >= self.identity and\
                            row['slen'].values[0] >= row['length'].values[0]:
                        # full length partial match
                        rank_list[scheme]['score'] += 70.0 / N
                        rank_list[scheme]['scheme'][locus] = \
                            '~{}'.format(row['number'].values[0])
                    elif row['coverage'].values[0] >= self.coverage and \
                            row['identity'].values[0] >= self.identity:
                        # partial length partial match
                        rank_list[scheme]['score'] += 20.0 / N
                        rank_list[scheme]['scheme'][locus] = \
                            '{}?'.format(row['number'].values[0])
                    else:
                        rank_list[scheme]['score'] += 0
                        rank_list[scheme]['scheme'][locus] = '-'
                else:
                    # multiples hits
                    for index, r in row.iterrows():
                        # if r['slen'] < r['length']:
                        #     if locus not in rank_list[scheme]['scheme']:
                        #         rank_list[scheme]['score'] += 20.0 / N
                        #         rank_list[scheme]['scheme'][locus] = \
                        #             '{}?'.format(r['number'])
                        #     else:
                        #         rank_list[scheme]['score'] -= 20.0 / N
                        #         rank_list[scheme]['score'] += 20.0 / N / len(row)
                        #         rank_list[scheme]['scheme'][locus] += \
                        #             '|' + '{}?'.format(r['number'])
                        if r['coverage'] == 1 and r['identity'] == 1:
                            # perfect match
                            if locus not in rank_list[scheme]['scheme']:
                                rank_list[scheme]['score'] += 100.0 / N
                                rank_list[scheme]['scheme'][locus] = \
                                    r['number']
                            else:
                                rank_list[scheme]['score'] -= 100.0 / N
                                rank_list[scheme]['score'] += 100.0 / N / len(row)
                                rank_list[scheme]['scheme'][locus] += \
                                    '|' + r['number']
                        elif r['coverage'] == 1 and\
                                r['identity'] >= self.identity and\
                                r['slen'] >= r['length']:
                            # full length partia match
                            if locus not in rank_list[scheme]['scheme']:
                                rank_list[scheme]['score'] += 70.0 / N
                                rank_list[scheme]['scheme'][locus] = \
                                    '~{}'.format(r['number'])
                            # else:
                            #     rank_list[scheme]['score'] -= 70.0 / N
                            #     rank_list[scheme]['score'] += 70.0 / N / len(row)
                            #     rank_list[scheme]['scheme'][locus] += \
                            #         '|' + '~{}'.format(r['number'])
                                # self.contamination = True
                        elif r['coverage'] >= self.coverage and\
                                r['identity'] >= self.identity:
                            # partial length partial match
                            if locus not in rank_list[scheme]['scheme']:
                                rank_list[scheme]['score'] += 20.0 / N
                                rank_list[scheme]['scheme'][locus] = \
                                    '{}?'.format(r['number'])
                            # else:
                                # rank_list[scheme]['score'] -= 20.0 / N
                                # rank_list[scheme]['score'] += 20.0 / N / len(row)
                                # rank_list[scheme]['scheme'][locus] += \
                                #     '|' + '{}?'.format(r['number'])
                                # self.contamination = True
                        else:
                            rank_list[scheme]['score'] += 0
                            if locus not in rank_list[scheme]['scheme']:
                                rank_list[scheme]['scheme'][locus] = '-'

        sorted_rank_list = sorted(rank_list.items(),
                                  key=lambda x: (x[1]['score']), reverse=True)
        bestscore = sorted_rank_list[0]
        self.scheme = bestscore[0]
        self.score = bestscore[1]


def main():
    print('Get out of here you little human!')

if __name__ == '__main__':
    main()
