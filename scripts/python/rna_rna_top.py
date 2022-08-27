import pandas as pd
import numpy as np
import sys
from collections import Counter

# RNA-DNA clusters
h5_clusters_dna = sys.argv[1]
# RNA only clusters
h5_clusters_rna = sys.argv[2]
# Out prefix
outprefix = sys.argv[3]

################################
# Functions
################################


def get_top_genes(h5_clusters_list, top_number, biotype, column_type):
    '''
    Count the frequency of RNA reads
    Args:
        h5_clusters_list(list):list of paths to h5 files (1 or more)
        top_number(int): Number of top expressing genes to return
        biotype(str): gene biotype (ie. protein_coding)
        column_type(str): "Exon" or "Intron"
    '''
    total_counts = Counter()
    for cluster_file in h5_clusters_list:
        total_counts.update(enumerate_columns_single(cluster_file, column_type))
    if biotype:
        filtered_counts = [(name, count) for (name, count) in total_counts.items() if biotype in name]
        total_counts = Counter(dict(filtered_counts))
    top_genes = total_counts.most_common(top_number) 
    return top_genes


def enumerate_columns_single(h5_clusters, column):
    '''
    Count elements in a column of a table
    Args:
        h5_clusters(str): Path to h5 file
        column(str): Column in table (e.g. 'Exon' or 'Intron')
    Returns:
        Counter
    '''
    if column == 'Exon':
        query = 'Exon != nan & Intron = nan & Repeat = nan'
    elif column == 'Intron':
        query = 'Exon = nan & Intron != nan & Repeat = nan'
    store = pd.HDFStore(h5_clusters)
    col = store.select('RPM', where=query, columns=[column])
    store.close()
    return Counter(list(col[column]))


#########################################################
# Function calls
#########################################################
h5_cluster_list = [h5_clusters_dna, h5_clusters_rna]
top_num = 10000
biotype = 'protein_coding'
column_type = 'Intron'

results = get_top_genes(h5_cluster_list, top_num, biotype, column_type)
df = pd.DataFrame(results, columns = ['Gene', 'Count'])
df.to_csv(outprefix + '.counts', sep='\t', header=True, index=False)

