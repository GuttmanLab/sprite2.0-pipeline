
import Levenshtein
from collections import defaultdict
import tqdm
from numba import njit, prange

seq_dict = defaultdict()
with open('/mnt/data/mm10_repeats_consensus.fasta') as in_fa:
    for name, seq in fasta_parse(in_fa):
        seq_dict[name] = seq


seqs = list(seq_dict.values())
len(seqs)

@njit(parallel=True)
def levenshtein_adj_list(seqs):
    ''' identify all umis within the hamming distance threshold (1 is best)
    and where the counts of the first umi is > (2 * second umi counts)-1
    will have duplicates'''

    adj_list = list()
    for i in prange(len(seqs)):
        dist_values = list()
        for j in range(i+1,len(seqs[i])):
            seq2 = seqs[j] #dict_keys object doesn't support indexing
            dist_values.append(Levenshtein.distance(seq, seq2))

        adj_list.append(dist_values)
    return adj_list


sparse_dict = levenshtein_adj_list(seqs)



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
