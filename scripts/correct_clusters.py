#%%
from collections import defaultdict
import re
import tqdm


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
    yield name, seq    

#%%
repeats = '/groups/guttman/Peter/genomes/GRCm38.p6/repeats/repetitive_sequences_20190828.fa'

correct_names = defaultdict(str)
with open(repeats) as in_fa:
    for name, seq in fasta_parse(in_fa):
        if 'tRNA' in name:
            correct_names[name.strip('>')] = name.strip('>').split(',')[0]


correct_names['Terra(14repeats)'] = 'Terra_14repeats'

#%%
cluster_file = "/groups/guttman/SPRITE/2019_09_12_NovaSeq_MegaSPRITE_Peter/HL522DSXX.all.clusters"
cluster_out = "/groups/guttman/SPRITE/2019_09_12_NovaSeq_MegaSPRITE_Peter.HL522DSXX.all.corrected.clusters"
with open(cluster_file) as in_cluster, \
open(cluster_out, 'w') as out_cluster:
    for cluster in tqdm.tqdm(in_cluster):
        # barcode, *reads = cluster.split('\t')
        corrected_cluster = cluster
        if 'tRNA' in cluster or 'Terra' in cluster:
            for k, v in correct_names.items():
                if k in cluster:
                    corrected_cluster = corrected_cluster.replace(k, v)
        out_cluster.write(corrected_cluster)
        


#%%
