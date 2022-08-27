import pandas as pd
import numpy as np
import sys
from collections import defaultdict
from itertools import combinations
import tqdm
import cooltools

# Clusters with DNA and RNA
h5_clusters_dna = sys.argv[1]
# Clusters with RNA Only
h5_clusters_rna = sys.argv[2]
# Bed File with Genes
filename = sys.argv[3]
# Prefix for saving contact matrix
outprefix = sys.argv[4]

#############################
# Scripts
#############################

def query_rnas_iteratively(h5_clusters, selection, column, chunksize):
    '''
    Read in chunks of h5 file query and process
    Args:
        h5_clusters(str): Path to h5 cluster file
        selection(list): gene names to search for
        column(str): Exon or Intron
        chunksize(int): Size of chunks to read one at a time
    '''
    store = pd.HDFStore(h5_clusters)
    if column == 'Exon':
        query = 'Exon=selection & Intron=nan & Repeat=nan'
    elif column == 'Intron':
        query = 'Exon=nan & Intron=selection & Repeat=nan'
    query_iter = store.select('RPM', where=query, columns=['Name', column], chunksize=chunksize)
    interactions = defaultdict(set)
    for chunk in tqdm.tqdm(query_iter):
        gb = chunk.groupby('Name')
        for g_name, g_df in gb:
            interactions[g_name].update(set(g_df[column]))
    store.close()
    return interactions


def make_matrix(interactions, genes):
    '''
    Make RNA-RNA contact matrix from dictionary of interactions
    Args:
        interactions(dict): interaction dictionary [barcode -> set(gene names)]
        genes(df): bedlike file of genes with column 'Gene' in correct order
    '''
    heatmap = defaultdict(lambda: defaultdict(int))
    for barcode, reads in tqdm.tqdm(interactions.items()):
         for r1,r2 in combinations(reads, 2):
             heatmap[r1][r2] += 1
             heatmap[r2][r1] += 1
    column_order = list(genes['Gene'])
    for gene in column_order: #Set self-interactions to zero
        heatmap[gene][gene] = 0
    df = pd.DataFrame.from_dict(heatmap)
    df = df[column_order].T
    df = df[column_order]
    return df

def ice_symmetric_matrix(df):
    '''
    ICE normalize an symmetric matrix
    Args:
       df(dataframe, np-array): matrix to be ICED
    '''
    normed  = cooltools.numutils.iterative_correction_symmetric(np.array(df))
    print(normed[2])
    new_df = pd.DataFrame(normed[0], index=df.index, columns =df.columns)
    return new_df

###############################
#Function Calls
###############################
chunk_size = 1000000
column_type = 'Intron'
genes = pd.read_csv(filename, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Gene'])
gene_list = list(genes['Gene'])

#Get RNA-RNA interactions in RNA-DNA clusters
interactions_dna = query_rnas_iteratively(h5_clusters_dna, gene_list, column_type, chunk_size)
matrix_dna = make_matrix(interactions_dna, genes) 
#Get RNA-RNA interactions in RNA Only clusters
interactions_rna = query_rnas_iteratively(h5_clusters_rna, gene_list, column_type, chunk_size)
matrix_rna = make_matrix(interactions_rna, genes)
#Merge matrices
matrix_total = matrix_dna.fillna(0)+matrix_rna.fillna(0)
#Ice Normalize
matrix_iced = ice_symmetric_matrix(matrix_total)

#Save Files
matrix_total.to_csv(outprefix + '.matrix', sep='\t', header=True, index=True)
matrix_iced.to_csv(outprefix + 'iced.matrix', sep='\t', header=True, index=True)

