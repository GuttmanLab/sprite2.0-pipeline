import pandas as pd
import numpy as np
import sys
from collections import defaultdict
from tqdm import tqdm
from itertools import combinations

# Path to DNA-RNA clusters
h5_clusters = sys.argv[1]
# Genes to select clusters
gene_file = sys.argv[2]
# Chromosome Sizes 
chrom_sizes_file = sys.argv[3]
# Outprefix to save final matrix
outprefix = sys.argv[4]


#########################################################
#Functions
#########################################################


def query_selected_dnas(h5_clusters, barcodes, chunksize):
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


def get_selected_weights(h5_clusters, barcodes):
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

def make_matrix(heatmap, interactions, weights):
    '''
    Make contact matrix, weighted by n/2
    Args:
        heatmap: defaultdict(lambda: defaultdict(int))
        interactions: Cluster dictionary(dict) = barcode -> set(reads)
        weights: Cluster sizes(dict) = barcode -> int(size)
    '''    
    for barcode, reads in tqdm(interactions.items()):
         count = 2/weights[barcode]
         for r1,r2 in combinations(reads, 2):
             heatmap[r1][r2] += count
             heatmap[r2][r1] += count
    return heatmap


def get_intron_clusters(h5_clusters, selection):
    '''
    Get barcodes for clusters containing intron reads of selected genes
    Args:
        h5_clusters(str): Path to RPM-DPM clusters
        selection(list): List of gene names
    '''
    query = 'Exon=nan & Intron=selection & Repeat=nan'
    store = pd.HDFStore(h5_clusters)
    col = store.select('RPM', where=query, columns=['Name', 'Intron'])
    store.close()
    return list(col['Name'])

def fill_matrix_labels(heatmap, labels):
    '''
    Fill in missing interactions from a interaction dictionary
    Args: 
        heatmap(dict): dictionary of interactions
        labels(list): list of interacting regions
    '''
    for r1,r2 in combinations(labels, 2):
         a = heatmap[r1][r2]
         b = heatmap[r2][r1]
    return heatmap

########################################################
#Function calls
########################################################
resolution = 1000000
chunk_size = 1000000

chrom_sizes = pd.read_csv(chrom_sizes_file, sep='\t', header=None)
genes = pd.read_csv(gene_file, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Gene'])
gene_list = list(genes['Gene'])

#Make ordered labels
labels = []
for index, row in chrom_sizes.iterrows():
    chrom = row[0]
    length = int(row[1])
    chrom_bins = [chrom + "_" + str(int(binned)) for binned in np.arange(0,length/resolution)]
    labels.extend(chrom_bins)

#Get clusters and reads
intron_barcodes = get_intron_clusters(h5_clusters, list(genes['Gene']))
intron_dna_clusters = query_selected_dnas(h5_clusters, intron_barcodes, chunk_size)
intron_dna_clusters_binned = bin_dna(intron_dna_clusters, resolution/1000)
intron_weights = get_selected_weights(h5_clusters, intron_barcodes)

#Make heatmap
heatmap = defaultdict(lambda: defaultdict(int))
heatmap = make_matrix(heatmap, intron_dna_clusters_binned, intron_weights)
heatmap_filled = fill_matrix_labels(heatmap, labels)
df = pd.DataFrame.from_dict(heatmap_filled)
df = df[labels]
df = df.T[labels]
print(df.sum().sum())

#Save matrix
df.to_csv(outprefix + '.matrix', sep='\t', header=True, index=True)

