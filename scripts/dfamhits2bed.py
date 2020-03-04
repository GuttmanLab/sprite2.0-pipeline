#%%
import gzip
from collections import defaultdict
import pybedtools
import subprocess
import re

#%%
anno_path = '/mnt/data/genomes/GRCm38.p6/repeats/mm10_dfam.nrph.hits.gz'
ref_fasta = '/mnt/data/genomes/GRCm38.p6/Mus_musculus.GRCm38.dna.primary_assembly.fa'

repeats_bed = dfam2bed(anno_path)


#%%
class bed:
    '''BED records holder

    Functions:
        add: Add new bed records to list
        tostring: convert entire list of records to string
    '''
    def __init__(self):
        self.chrom = []
        self.start = []
        self.end = []
        self.strand = []

    def add(self, chrom, start, end, strand=None):
        self.chrom.append(chrom)
        self.start.append(start)
        self.end.append(end)
        if strand != None:
            self.strand.append(strand)

    def tostring(self, chr_format=True):
        out_str = ''

        for i in range(0, len(self.chrom)):
            if chr_format:
                out_str += '>' + str(self.chrom[i]) + '\t' + str(self.start[i]) + \
                '\t' + str(self.end[i]) + '\n'
            else:
                out_str += str(self.chrom[i].split('chr')[-1]) + '\t' + str(self.start[i]) + \
                '\t' + str(self.end[i]) + '\n'

        return out_str

#%%
test = bed()
test.add('chr1', 100, 234, '+')


#%%

def dfam2bed(anno_path):
    '''Convert dfam non-redundent hits annotation into a bed file, order by family 
    and print stats for each family

    Args:
        anno_path(str): path to dfam nrph annotation file
        bed_out(str): path of output bed file

    Notes:

        Dfam.hits
        =========
        Dfam.hits file contains a tab separated list of the matches in dfamseq to each 
        of the models contained in Dfam.hmm that score above the GA threshold for the 
        model.  This files contains all hits found by the Dfam HMMs and does not try 
        and resolve overlaps between related models. The fields in that file are as 
        follows:

        Dfamseq identifier   - Constructed from organism name, chromosome (or NC for 
                            unplaced contigs), sequence number.
        Dfam accession       - i.e. DFXXXXXXX, use this to link on.
        Dfam identifier      - 16 character name for the entry.
        bit score            - The score of the match
        E-value              - Estimate of the statistical significance of this match in 
                            context of the the size of dfamseq.
        bias                 - Bias composition of the hit
        model start          - The first position matched in the model by the hit.
        model end            - The last position matched in the model by the hit.
        strand               - Either '+' or '-'.
        alignment start      - The position of the first base in the alignment between 
                            the hit and the model.
        alignment end        - The start position of the first base in the alignment 
                            between the hit and the model.
        envelope start       - The start position of the sequence that defines a 
                            subsequence for which there is substantial probability mass 
                            supporting a homologous hit.  See the HMMER documentation for
                            further information on HMM envelopes.
        envelope end         - The end position of the sequence that defines a subsequence 
                            for which there is substantial probability mass supporting 
                            a homologous hit.
        sequence length      - The length of the sequence containing the match.
        Kimura divergence    - The Kimura divergence of the aligned bases

    '''
    #order by family
    fam_dict = defaultdict(bed)
    # count = 0
    with gzip.open(anno_path, 'r') as dfam_anno:
        header = dfam_anno.readline()
        for line in dfam_anno:
            chrom, _, family, _, _, _, _, _, strand, start, end, *rest = \
            line.decode('utf-8').rstrip('\n').split('\t')
            # count += 1
            if strand == '-':
                fam_dict[family].add(chrom, end, start, strand)
            elif strand == '+':
                fam_dict[family].add(chrom, start, end, strand)
            else:
                fam_dict[family].add(chrom, start, end)

    return fam_dict    



#%%
def fasta_parse(fp):
    '''Parse fasta files

    Args:
        fp (str): Open fasta file
    '''

    record_out = 0
    record_count = 0
    name, seq = [''] * 2

    for line in fp:

        try:
            rec = line.decode('utf-8').rstrip()
        except AttributeError:
            rec = line.rstrip()

        if rec.startswith('>'):
            if record_count - record_out > 0:
                 yield name, seq
                 name, seq = [''] * 2

            name = rec
            record_count += 1
        elif line.startswith('@'):
            raise Exception('Input starts with @ suggesting this is a fastq no fasta')
        else:
            seq += rec


#%%
def get_sequence(bed, fasta):
    '''Get sequences of coordinates in bed file
    pybedtools, get sequences from genome file for each family

    Args:
        bed(class): bed class 
        fasta(str): path to fasta file
    '''
    



#make consensus with degenerate bases for families


#%%
a = pybedtools.BedTool("""
chr1 1 10
chr1 50 55""", from_string=True)

#%%
# B1_Mm
# 'MER70C'
#loop through repeats
for k, v in repeats_bed:
    bed_file = v.tostring()

bt_bed = pybedtools.BedTool(repeats_bed['B1_Mm'].tostring(chr_format=False), 
from_string=True)
bt_seqs = bt_bed.sequence(fi=ref_fasta)
seqs = open(bt_seqs.seqfn).read()

fasta_str = ''
for name, seq in fasta_parse(open(bt_seqs.seqfn)):
        fasta_str += name + '\n' + seq + '\n'
    
#%%
kalign_msa(fasta_str)

# '/mnt/data/mer70c.fasta'
with open('/mnt/data/B1_Mm.fasta', 'w') as out_fa:
    out_fa.write(fasta_str)

#%%
def kalign_msa(fasta):
    '''Multiple sequence alignment for read loss analysis
    :param dict seq_counter_dict: dict of umi:Counter(sequences) object with sequences
    :param qual_dict: perserve fq quality of aligned reads
    '''
    
    #Multiple sequence alignment
    #http://msa.sbc.su.se/cgi-bin/msa.cgi
    kalign_cmd = ['kalign', '-f', 'fasta'] #'-s', '80.0', '-e', '3.0', '-t', '3.0', '-m', '0.0'


    p = subprocess.Popen(kalign_cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
    kalign_stdout = p.communicate(input=fasta_str.encode('utf-8'))[0]

    #parse kalign output and derive consensus seq
    head, sep, tail = kalign_stdout.partition(b'>') #remove head kalign intro

    alignment = sep.decode() + tail.decode()
    alignment = re.sub('\n', '', alignment)
    alignment = list(filter(None, re.split('(>\d+:\d+-\d+)', alignment)))

    with open('/mnt/data/mer70c.txt', 'w') as out:
        for qname, seq in fasta_parse(alignment):
            out.write(seq + '\n')


    return alignment




#%%
