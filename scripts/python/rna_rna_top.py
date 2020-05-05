import pandas as pd
import numpy as np
from collections import Counter
from timeit import Timer
import matplotlib.pyplot as plt


def main():
    h5_clusters = '/mnt/data/RNA_DNA_SPRITE/20200315_RNA_DNA_Mario/8226.comboall.blosc.h5'
    # store = pd.HDFStore(h5_clusters)
    # store.groups()
    # store.keys()


    select = get_top_genes(h5_clusters, 50)

    t1 = Timer(lambda: enumerate_columns_whole(h5_clusters, 'exon'))
    print(t1.timeit(number=1)) #4.892679850017885

    t2 = Timer(lambda: enumerate_columns_exon(h5_clusters))
    print(t2.timeit(number=1)) #25.765680634009186




def get_top_genes(h5_clusters, count_filter):
    '''
    Find most frequent genes for RNA-RNA plotting

    Args:
        h5_clusters(str): Path to h5 file
        count_filter(int): Cutoff for gene count
    '''

    exon_counts = enumerate_columns_whole(h5_clusters, 'exon')
    intron_counts = enumerate_columns_whole(h5_clusters, 'intron')

    combined_counts = exon_counts + intron_counts
    # plt.hist(combined_counts.values(), bins=[1, 10, 20, 30, 40, 50, 100, 500])

    top_genes = [k for k, v in combined_counts.items() if v > count_filter]
    #remove intron exon distinction 
    top_out = set([i.split('.')[0] for i in top_genes])

    return top_out



# def enumerate_introns():
#     '''
#     Get count of genes in clusters
#     '''

#     interactions = defaultdict(set)

#     for df in pd.read_hdf(h5_clusters, 'RPM', where='exon != NA & intron != NA)',
#                      columns=['Name', 'exon', 'intron'], chunksize=1000000)
        
#     return 

# if with_dna:
#         #only select RNA from clusters also containing DNA
#         dpm_clusters = column_groups(h5_clusters, 'DPM', 'Name')


def enumerate_columns_whole(h5_clusters, column):
    '''
    Count elements in a column of a table

    Args:
        h5_clusters(str): Path to h5 file
        column(str): Column in table (e.g. exon)

    Returns:
        Counter
    '''

    store = pd.HDFStore(h5_clusters)
    col = store.select_column('RPM', column)
    store.close()
    values = col[np.logical_not(pd.isnull(col))]
    groups = Counter(list(values))
    return groups


def enumerate_columns_exon(h5_clusters):
    '''
    Count elements in a column of a table

    Args:
        h5_clusters(str): Path to h5 file
        column(str): Column in table (e.g. exon)

    Returns:
        Counter
    '''

    col = pd.read_hdf(h5_clusters, 'RPM', where='exon != NA', columns=['exon'])
    groups = Counter(list(col))
    
    return groups



def column_groups(h5_clusters, table, column):
    '''
    Get unique (groups) of a column in a table

    Args:
        h5_clusters(str): Path to h5 file
        table(str): Either RPM or DPM table
        column(str): Column in table (e.g. Name)
    '''
    store = pd.HDFStore(h5_clusters)
    groups = store.select_column(table, column).unique()
    store.close()    
    return groups
