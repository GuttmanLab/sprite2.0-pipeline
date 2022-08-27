import pandas as pd
import numpy as np
import sys
from tqdm import tqdm
from collections import defaultdict


# Path to DNA-RNA clusters
h5_clusters = sys.argv[1]
# Genes to select clusters
gene_file = sys.argv[2]
# Outprefix to save final bedgraphs
outprefix = sys.argv[3]

##########################################################################
## Functions
##########################################################################


def query_rnas_iteratively(h5_clusters, selection, chunksize):
    '''
    Read in chunks of h5 file query and process
    Args:
        h5_clusters(str): Path to h5 cluster file
        selection(list): gene names to search for
        chunksize(int): Size of chunks to read in one at a time
    '''
    store = pd.HDFStore(h5_clusters)
    query = 'Intron=selection & Exon=nan & Repeat=nan'
    query_iter = store.select('RPM', where=query, columns=['Name', 'Intron'], chunksize=chunksize)
    interactions = defaultdict(set)
    for chunk in tqdm(query_iter):
        gb = chunk.groupby('Name')
        for g_name, g_df in gb:
            interactions[g_name].update(set(g_df['Intron']))
    store.close()
    return interactions

def query_dnas_iteratively(h5_clusters, barcodes, chunksize):
    '''
    Get DPM reads for a selection of barcodes. Deduplicate DPM reads by binning to 1000bp
    Args:
        h5_clusters(str): Path to RPM-DPM clusters
        barcodes(list): List of cluster barcodes to query
        chuncksize(int): size of chunks to read
    Returns:
        Cluster dictionary(dict): barcode -> set(reads), read format = 'chrX_100'
    '''
    interactions = defaultdict(set)
    store = pd.HDFStore(h5_clusters)
    query_iter = store.select('DPM', where='Name=barcodes', columns=['Name', 'Chromosome', 'Start'], chunksize=chunksize)
    for chunk in tqdm(query_iter):
        chunk['Key'] = chunk['Chromosome'] + '_' + chunk['Start'].floordiv(1000).astype(str)
        mini_dict = chunk.groupby('Name', sort=False)['Key'].apply(set).to_dict()
        interactions.update(mini_dict)
    store.close()
    return interactions

def bin_dna(interactions, resolution):
    '''
    Bin DNA to a given resolution, keeping multiple reads per bin
    Args:
         interactions(dict): Cluster dictionary {barcode -> set(reads binned at 1000bp)}
         resolution(int): Desired resolution/1000
    Returns:
         Cluster dictionary(dict): barcode -> list (reads binned at resolution)
    '''
    interactions_binned = defaultdict(list)
    for barcode, reads in tqdm(interactions.items()):
        binned = [read.split('_')[0] + '_' + str(int(int(read.split('_')[1])//resolution)) for read in reads]
        interactions_binned[barcode] = binned
    return interactions_binned

def get_cluster_weights(h5_clusters, barcodes):
    '''
    Get cluster sizes for a selection of clusters
    Args:
         h5_clusters(str): Path to RPM-DPM clusters
         barcodes(list): list of barcodes to query
    '''
    store = pd.HDFStore(h5_clusters)
    query = store.select('Names', where='Name=barcodes', columns = ['Name', 'Size'])
    store.close()
    return query.set_index('Name')['Size'].to_dict()


def make_bedgraph(dna_clusters, rna_clusters, weights):
    '''
    Make RNA-DNA interaction matrix
    Args:
        dna_clusters(dict):barcode -> binned DNA reads
        rna_clusters(dict):barcode -> RNA reads
        weights(dict):barcode -> cluster size
    '''
    bed = defaultdict(lambda: defaultdict(int))
    for key, value in dna_clusters.items():
        rna_names = rna_clusters[key]
        count = 2/weights[key]
        for rna in rna_names:
            for read in value:
                bed[rna][read]+=count
    df = pd.DataFrame.from_dict(bed)
    return df


########################################################
#Function calls
########################################################
resolution = 1000000
chunk_size = 1000000

genes = pd.read_csv(gene_file, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Gene'])
gene_list = list(genes['Gene'])


interactions_rna = query_rnas_iteratively(h5_clusters, gene_list, chunk_size)
barcodes = list(interactions_rna.keys())
interactions_dna = query_dnas_iteratively(h5_clusters, barcodes, chunk_size)
interactions_dna_binned = bin_dna(interactions_dna, resolution/1000)
weights = get_cluster_weights(h5_clusters, barcodes)

bedgraphs = make_bedgraph(interactions_dna_binned, interactions_rna, weights)

bedgraphs.to_csv(outprefix + '.matrix', sep='\t', header=True, index=True)

