import pandas as pd
from collections import Counter, defaultdict
import pyranges as pr
from scipy.sparse import coo_matrix, dok_matrix, save_npz, load_npz
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
import gc
import json
import os
import argparse
from rna_rna_top import get_top_genes


from timeit import Timer
from memory_profiler import profile
from cluster import file_open


def parse_arguments():
    parser = argparse.ArgumentParser(description =
            "Clusters to other format conversion")
    parser.add_argument('-i', '--input', action = 'store', metavar = 'FILE',
                        help = 'h5 file')
    parser.add_argument('-o', '--output', action = 'store', metavar = 'DIRECTORY',
                        help = 'Output directory') 
    parser.add_argument('--gtf', metavar = 'GTF', action = 'store',
                        help = "GTF from which to retrieve gene names")
    parser.add_argument('--biotype', metavar = 'BIOTYPE', action = 'store',
                        default = 'protein_coding', 
                        help = "Biotype of genes") 
    parser.add_argument('--region', metavar = 'REGION', type = str,
                        action = 'store', default = 'exons',
                        choices=['exons', 'introns'],
                        help = "Interactions of either exons or introns")
    parser.add_argument('--chromosome', metavar = 'CHROM', action = 'store',
                        default = None, 
                        help = 'Get interactions of genes on only a single \
                                chromosome or genome wide') 
    parser.add_argument('--chunksize', metavar = 'CHUNK', action = 'store',
                        default = 1000000, 
                        help = "Size of chunks to be read in one at a time") 
    parser.add_argument('--self', action = 'store_true', 
                        help = "Include self interactions (singletons)")
    parser.add_argument('--top', metavar = 'INT', action = 'store', type = int,
                    default = 0, 
                    help = "Count cutoff of genes in clusters to be used for rna-rna \
                            interaction")
                       
    return parser.parse_args()


def main():

    args = parse_arguments()

    gene_data = genes(args.gtf, args.output, args.region, args.top)
    gene_data.get_gene_names(chrom=args.chromosome, biotype=args.biotype)
    gene_data.write_to_file()

    if args.top > 0:
        top_genes = get_top_genes(args.input, args.top) #returns a set
        gene_data._gene_names = [i for i in gene_data._gene_names if i in top_genes]
        gene_data.write_to_file(top=True)

    if args.region == 'introns':
        select_introns = gene_data.to_introns()
        interactions = query_introns_iteratively_read(args.input, select_introns, 
                                                      args.chunksize)
        gene_data.write_interactions(interactions)
        gene_data.make_matrix(interactions, 'introns', args.self, False)

    elif args.region == 'exons':
        select_exons = gene_data.to_exons()
        interactions = query_exons_iteratively_read(args.input, select_exons, 
                                                    args.chunksize)
        gene_data.write_interactions(interactions)
        gene_data.make_matrix(interactions, 'exons', args.self, False)

    # S_coo = load_npz('/mnt/data/RNA_DNA_SPRITE/20200315_RNA_DNA_Mario/gene_lists/2_protein_coding.npz')


class genes():
    '''
    Holder for gene names used for searches of HDF5 file
    '''

    def __init__(self, gtf, output_dir, region, top=0):
        self._gtf = gtf
        self._gene_names = []
        self._output_dir = output_dir
        self._region = region
        self._top = top

    def get_gene_names(self, chrom=None, biotype='protein_coding'):
        '''
        Get protein coding (or other biogroup) genome wide or specific chromosome

        Args:
            self.gft(str): Path to gtf from which to extract gene names
            chrom(str): Which chromosome to extract gene names from (ensembl names)
            biotype(str): Biotype of gene to extract
        '''
        if chrom == None:
            self._chrom = 'All'
        else:
            self._chrom = chrom

        self._biotype = biotype

        #Check if output file already exists, and if it does, read from it
        bn = os.path.basename(self._gtf).split('.gtf.gz')[0]
        out_file = self._output_dir + '/' + bn + '_' + self._chrom + '_' + self._biotype + '.txt'

        if os.path.isfile(out_file):
            self.read_from_file(out_file)
        else:

            gtf_gr = pr.read_gtf(self._gtf)
            gtf_gr = gtf_gr[gtf_gr.Feature == 'gene']
            gtf_gr = gtf_gr[gtf_gr.gene_biotype == biotype]
            if chrom != None:
                #Get genes on chromosome
                gtf_gr = gtf_gr[gtf_gr.Chromosome == str(chrom)]
            #ensure genome order (pyranges seems to always sort by strand)
            gtf_df = gtf_gr.df.sort_values(by=['Chromosome', 'Start'])

            unique_gene_name = unique(list(gtf_df['gene_name']))
            self._gene_names = unique_gene_name
            
            self._gtf_gr = gtf_gr


    def to_introns(self):
        return [i + '.intron' for i in self._gene_names]


    def to_exons(self):
        return [i + '.exon' for i in self._gene_names]

    def to_index(self, g_type):
        '''
        Get exon or intron index used for making matrix entries 

        Args:
            g_type(str): Make index for introns or exons
        '''
        if g_type == 'introns':
            t = '.intron'
        elif g_type == 'exons':
            t = '.exon'

        idx = defaultdict(int)
        j = 0
        for i in self._gene_names:
            idx[i + t] = j
            j += 1

        return idx

    def write_to_file(self, top=False):
        '''
        Write gene_names to file
        '''
        bn = os.path.basename(self._gtf).split('.gtf.gz')[0]
        if top:
            out_file = self._output_dir + '/' + bn + '_' + self._chrom + '_' + \
                       self._biotype + '_top' + str(self._top) + '.txt'
        else:
            out_file = self._output_dir + '/' + bn + '_' + self._chrom + '_' + \
                       self._biotype + '.txt'
            
        if os.path.isfile(out_file):
            print('File already exists')
        else:
            with open(out_file, 'w') as out:
                for i in self._gene_names:
                    out.write(i + '\n')

    def read_from_file(self, in_file):
        '''
        Read gene_names from file
        '''
        
        with open(in_file) as read:
            for i in read:
                self._gene_names.append(i.rstrip())

    def write_interactions(self, interactions):
        '''
        Write out interactions dictionary
        '''
        if self._top > 0 :
            out_file = self._output_dir + '/' + self._chrom + '_' + self._biotype + '_' + \
                   self._region + '_top' + str(self._top) +'.npy'
        else:
            out_file = self._output_dir + '/' + self._chrom + '_' + self._biotype + '_' + \
                    self._region + '.npy'
        np.save(out_file, interactions)


    def make_matrix(self, interactions, g_type, self_interaction=False, plot=False):
        '''
        Make matrix from dictionary

        Args:
            interactions(dict): gene interaction dictionary 
            self_interactions(logical): ignore self or include
        '''

        gene_idx = self.to_index(g_type)

        S = dok_matrix((len(gene_idx.keys()),len(gene_idx.keys())), dtype='float32')
        count_out = 0
        for barcode, items in interactions.items():
            if not self_interaction:
                if len(items) > 1:
                    idx = [gene_idx.get(i) for i in items]
                    for i, j in combinations(idx, 2):
                        S[i,j] += 1
                        S[j,i] += 1
                        count_out += 1
            elif self_interaction:
                if len(items) <= 1:
                    idx = gene_idx.get(list(items)[0])
                    try:
                        S[idx, idx] += 1
                    except IndexError:
                        print('Index error:', idx)
                else:
                    idx = [gene_idx.get(i) for i in items]
                    for i, j in combinations(idx, 2):
                        S[i,j] += 1
                        S[j,i] += 1
                        count_out += 1
                    #add self
                    for x in idx:
                        S[x,x] += 1

        S_coo = S.tocoo()

        if self._top > 0 :
            out_file = self._output_dir + '/' + self._chrom + '_' + self._biotype + '_' + \
                   self._region + '_top' + str(self._top) +'.npz'
        else:
            out_file = self._output_dir + '/' + self._chrom + '_' + self._biotype + '_' + \
                    self._region + '.npz'

        save_npz(out_file, S_coo)

        if plot:
            S_dense = S_coo.todense()
            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(111)
            im = ax.matshow(np.log10(S_dense), cmap='YlOrRd')
            fig.colorbar(im)





    # def get_order(self):
    #     '''
    #     Get order of genes on chromosome

    #     Args:
    #         type(str): Either 'exons' or 'introns'
    #     '''      
        
    #     items = list(self.gene_names)
    #     items_idx = sorted(range(len(items)), key=lambda pos: self.index.get(items[pos]))
    #     self.order = items_idx
        


def unique(items):
    '''
    Deduplicate list, preserving order

    Args:
        items(list): List of items to deduplicate
    '''
    seen = set()
    seen_add = seen.add
    unq = [i for i in items if not (i in seen or seen_add(i))]
    return unq


def base_introns(bed_file, selection):
    '''
    Get introns from the raw bed format to assess benefit of h5 format

    Args:
        bed_file(str): Path to bed file
        selection(list): Which genes to select
    '''
    interactions = defaultdict(set)
    with file_open(bed_file) as in_file:
        for line in in_file:
            chrom, start, end, strand, barcode, score, exon, intron, repeat, r_type = \
                line.decode('utf-8').rstrip().split('\t')
            if r_type == 'RPM':
                if intron in selection and exon == 'NA':
                    interactions[barcode].update(intron)

    return interactions

# @profile
def query_introns(h5_clusters, selection):
    '''
    Read in entire h5 file query and process

    Args:
        h5_clusters(str): Path to h5 indexed file
        selection(list): Which genes to query
    '''

    query_df = pd.read_hdf(h5_clusters, 'RPM', 
                           where='exon = NA | exon = nan & intron = selection',
                           columns=['Name','intron'])

    #Iterate over groups(clusters) and create intron matrix
    gb = query_df.groupby('Name')

    interactions = defaultdict(set)
    for g_name, g_df in gb:
            interactions[g_name].update(set(g_df.intron))

    return interactions

#Create a dictionary of intron combinations by barcode
# @profile
def query_introns_select(h5_clusters, selection):
    '''
    Read in chunks of h5 file query and process

    Args:
        h5_clusters(str): Path to h5 indexed file
        selection(list): Which genes to query
        chunksize(int): Size of chunks to read in one at a time
    '''
    store = pd.HDFStore(h5_clusters)

    query_df = store.select('RPM', where='exon = NA | exon = nan & intron = selection', 
                 columns=['Name','intron'])

    #Iterate over groups(clusters) and create intron matrix
    gb = query_df.groupby('Name')

    interactions = defaultdict(set)
    for g_name, g_df in gb:
            interactions[g_name].update(set(g_df.intron))
    
    store.close()
    return interactions


# @profile
def query_introns_iteratively(h5_clusters, selection, chunksize):
    '''
    Read in chunks of h5 file query and process

    Args:
        h5_clusters(str): Path to h5 indexed file
        selection(list): Which genes to query
        chunksize(int): Size of chunks to read in one at a time
    '''
    store = pd.HDFStore(h5_clusters)

    query_iter = store.select('RPM', where='exon = NA | exon = nan & intron = selection', 
                 columns=['Name','intron'], chunksize=chunksize)

    interactions = defaultdict(set)
    for chunk in query_iter:
        gb = chunk.groupby('Name')
        for g_name, g_df in gb:
                interactions[g_name].update(set(g_df.intron))
    store.close()
    return interactions


def query_introns_iteratively_read(h5_clusters, selection, chunksize):
    '''
    Read in chunks of h5 file query and process

    Args:
        h5_clusters(str): Path to h5 indexed file
        selection(list): Which genes to query
        chunksize(int): Size of chunks to read in one at a time
    '''
    interactions = defaultdict(set)

    for chunk in pd.read_hdf(h5_clusters, 'RPM', where='exon = NA | exon = nan & intron = selection', \
                             columns=['Name','intron'], chunksize=chunksize):
        gb = chunk.groupby('Name')
        for g_name, g_df in gb:
                interactions[g_name].update(set(g_df.intron))
        
    return interactions


def query_exons_iteratively_read(h5_clusters, selection, chunksize):
    '''
    Read in chunks of h5 file query and process

    Args:
        h5_clusters(str): Path to h5 indexed file
        selection(list): Which genes to query
        chunksize(int): Size of chunks to read in one at a time
    '''
    interactions = defaultdict(set)

    for chunk in pd.read_hdf(h5_clusters, 'RPM', where='exon = selection', \
                             columns=['Name','exon'], chunksize=chunksize):
        gb = chunk.groupby('Name')
        for g_name, g_df in gb:
                interactions[g_name].update(set(g_df.exon))
        
    return interactions

def speed_tests(h5_clusters, select_introns):

    '''
    bed file parse
    72.89146137200441
    bed.gz file parse
    86.70955644599599

    lzo compression
    query_intron
    40.68881972900272
    query_intron_select
    40.34403066700179
    query_intron_iteratively_500000
    16.294188640000357
    query_intron_iteratively_100000
    16.731297678001283
    query_intron_iteratively_1000000
    16.45854562599561
    query_intron_iteratively_read_1000000
    16.37871296799858

    blosc compression
    query_intron
    40.76128340699506
    query_intron_select
    38.06138102099794
    query_intron_iteratively_500000
    12.413355795004463
    query_intron_iteratively_100000
    12.64127059099701
    query_intron_iteratively_1000000
    12.367498462001095
    query_intron_iteratively_read_1000000
    12.533995552003034

    No compression
    query_intron
    39.15648021600646
    query_intron_select
    38.52819720399566
    query_intron_iteratively_500000
    10.51807031899807
    query_intron_iteratively_100000
    10.285655829997268
    query_intron_iteratively_1000000
    9.788864490998094
    query_intron_iteratively_read_1000000
    9.666394640997169
    '''

    gtf = "/mnt/data/genomes/GRCh38/Homo_sapiens.GRCh38.97.gtf.gz"
    gene_data = genes(gtf, '/mnt/data/RNA_DNA_SPRITE/20200315_RNA_DNA_Mario/gene_lists')
    gene_data.get_gene_names(chrom='2')
    gene_data.write_to_file()

    h5_clusters = '/mnt/data/RNA_DNA_SPRITE/20200315_RNA_DNA_Mario/8226.comboall.blosc.h5'
    bed_clusters = '/mnt/data/RNA_DNA_SPRITE/20200315_RNA_DNA_Mario/8226.comboall.bed.gz'

    select_introns = gene_data.to_introns()
    select_exons = gene_data.to_exons()


    print('query_intron')
    t = Timer(lambda: query_introns(h5_clusters, select_introns))
    print(t.timeit(number=1))
    print('query_intron_select')
    t2 = Timer(lambda: query_introns_select(h5_clusters, select_introns))
    print(t2.timeit(number=1))
    print('query_intron_iteratively_500000')
    t3 = Timer(lambda: query_introns_iteratively(h5_clusters, select_introns, 500000))
    print(t3.timeit(number=1))
    print('query_intron_iteratively_100000')
    t4 = Timer(lambda: query_introns_iteratively(h5_clusters, select_introns, 100000))
    print(t4.timeit(number=1))
    print('query_intron_iteratively_1000000')
    t5 = Timer(lambda: query_introns_iteratively(h5_clusters, select_introns, 1000000))
    print(t5.timeit(number=1))
    print('query_intron_iteratively_read_1000000')
    t6 = Timer(lambda: query_introns_iteratively_read(h5_clusters, select_introns, 1000000))
    print(t6.timeit(number=1))
    t7 = Timer(lambda: base_introns(bed_clusters, select_introns))
    print(t7.timeit(number=1))

    

if __name__ == "__main__":
    main()
    # speed_tests()