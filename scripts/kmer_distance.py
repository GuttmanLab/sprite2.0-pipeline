#%%
import gzip
from collections import defaultdict
import functools
import numpy as np
import pandas as pd
import seaborn as sns

#%%
# generate_kmers('AGATGACGATAGACCAGATAGACAGATAG', 5)
repeats = '/mnt/data/genomes/GRCm38.p6/repeats/repetitive_sequences_20190828.fa'
min_len = shortest_seq(repeats)
print(min_len)
repeats_kmers = sequence2kmers(repeats, 50)
counts_dict = count_intersection(repeats_kmers)
counts_matrix = dict2matrix(counts_dict)

counts_matrix.to_excel('/mnt/data/genomes/GRCm38.p6/repeats/kmer_matrix_50mer.xlsx')

#%%
def count_intersection(kmers):
    '''Count number of shared kmers between sequences

    Args:
        kmers(dict): Dictionary of sequence names (keys) and kmers (values)

    Notes:
        s1 = {'ATGA', 'ATGC'}
        s2 = {'ATGA', 'ACCT'}

        s1.intersection(s2)
    '''
    count_dict = defaultdict(list)
    for k, v in kmers.items():
        for k_2, v_2 in kmers.items():
            shared_kmers = len(v.intersection(v_2))
            count_dict[k].append(shared_kmers)

    return count_dict

#%%
def dict2matrix(dictionary):
    '''Convert a counts dictionary to a numpy matrix for plotting

    Args:
        dictionary(dict): Counts matrix from count_intersection function
    '''
    matrix = np.array([v for v in dictionary.values()])
    
    names = list(dictionary.keys())

    out_matrix = pd.DataFrame(matrix, index=names, columns=names)

    return out_matrix

# dataDict = {'device1':(1,1,0,1), 'device2':(0,1,0,1), 'device3':(1,0,0,1)}
# orderedNames = ['device1','device2','device3']

# dataMatrix = np.array([dataDict[i] for i in orderedNames])



#%%
def shortest_seq(fasta):
    '''Gives length of shortes sequence in fasta
    Args:
        fasta(str): Path to fasta file
    '''
    lengths = []
    with file_open(fasta) as in_fa:
        for qname, seq in fasta_parse(in_fa):
            lengths.append(len(seq))

    return min(lengths)


#%%
def generate_kmers(sequence, size=50):
    '''Generate all possible kmers of size k from sequence
    
    Args:
        sequence(str): Sequence form which to create kmers
        size(int): Size of kmers to generate
    '''

    if len(sequence) < size:
        print('Sequence shorter than kmer size')
        return set(sequence)
    else:
        kmers = set()
        seq_len = len(sequence)
        position = 0
        while seq_len > size+position:
            kmers.add(sequence[position:size+position])
            position += 1
        return kmers

#%%
def expand_N(kmers):
    '''If N present in sequence, expand with all possible bases ATGC
        
    Args:
        kmers(dict): kmer dictionary with name (keys) and all possible kmers (values)    
    '''

    for k, v in kmers.items():



#%%
def sequence2kmers(fasta, kmer_size):
    '''From a fasta file generate a dictionary with all kmers for each fasta record

    Args:
        fasta(str): Path to fasta file
        kmer_size(int): size of kmers to generate from fasta records
    '''
    kmers = defaultdict()
    with file_open(fasta) as in_fa:
        for qname, seq in fasta_parse(in_fa):
            kmers[qname] = generate_kmers(seq, kmer_size)

    return kmers




#%%
def file_open(filename):
    '''
    Open as normal or as gzip
    Faster using zcat?

    Args:
        filename(str): Path to file to be opened
    '''
    #does file exist?
    f = open(filename,'rb')
    if (f.read(2) == b'\x1f\x8b'): #compressed alsways start with these two bytes
        f.seek(0) #return to start of file
        return gzip.GzipFile(fileobj=f, mode='rb')
    else:
        f.seek(0)
        return f

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
        elif rec.startswith('@'):
            raise Exception('Input starts with @ suggesting this is a fastq no fasta')
        else:
            seq += rec

#%%
