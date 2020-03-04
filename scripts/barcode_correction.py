
#cython: language_level=3

# import pyximport; pyximport.install()
# import _edit_distance
from collections import defaultdict, Counter
import re
import gzip
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm
# from scipy.spatial import distance
from babrahamlinkon.general import fastq_parse, file_open
from babrahamlinkon._dedup_umi import edit_distance
from copy import deepcopy
import argparse
import numba
from itertools import combinations

'''
1. convert barcodes into simple strings to allow hamming distance calculation
2. create dictionary to hold barcodes and corresponding qnames
3. collapse keys with 1 hamming distance
'''



# distance.hamming('ABC','BBC')


# edit_distance('ABC'.encode(), 'BBD'.encode())


#TODO:
# def read_config(path):




#TODO:
# class Barcodes():
#     '''Holder and manipulation of barcodes from BarcodeID
#
#     Expects a list of barcodes either in the long format:
#     [TermStag_bot_3][R1_bot_A1][R2_bot_B5][R1_bot_A4][R2_bot_B2][R1_bot_A9][R2_bot_B11][R1_A1_1000x]
#     or in the simple format:
#
#
#     '''
#
#     __init__(self, barcodes):
#
#         if '.' in barcodes[0]:
#             simple_barcode = TRUE
#         else:
#             simple_barcode = FALSE
#
#         self.barcodes = barcodes
#         self.simple_barcode = simple_barcode
#         self.simple_convert = {'1':'A', '2':'B', '3':'C', '4':'D', '5':'E', '6':'F', '7':'G', '8':'H', '9':'I',
#                                '10':'J', '11':'K', '12':'L', 'FOUND':'M',
#                                '1000x':'4', '100x':'3', '10x':'2', '1x':'1', 'empty':'0'}
#         self.long_convert =
#
#     def convert2simple(self):
#         '''Convert full barcodes to simple format
#         '''
#         if self.simple_barcode:
#             print('Barcodes already in simple format')
#             return
#         else:



# ['TermStag_bot_3',
#  'R1_bot_A1',
#  'R2_bot_B5',
#  'R1_bot_A4',
#  'R2_bot_B2',
#  'R1_bot_A9',
#  'R2_bot_B11',bc_list
#  'R1_A1_1000x']




# cluster_sizes = {'1000x':Counter(), '100x':Counter(), '10x':Counter(), '1x':Counter(), 'empty':Counter()}
# round = {1:Counter(),2:Counter(),3:Counter(),4:Counter(),5:Counter(),6:Counter(),7:Counter(),8:Counter()}


def get_barcodes(fastq):
    '''Get the barcodes from the qname of a fastq file

    Args:
        fastqgzfile (str): Path to the gzipped fastq file

    Return:
        A dictionary of simplified barcodes (key) and a list of their qnames (values)
    '''

    convert_dict = {'1':'A', '2':'B', '3':'C', '4':'D', '5':'E', '6':'F', '7':'G', '8':'H', '9':'I',
                    '10':'J', '11':'K', '12':'L', 'FOUND':'M',
                    '1000x':'4', '100x':'3', '10x':'2', '1x':'1', 'empty':'0'}

    skipped = 0
    total = 0
    bc_container = defaultdict(list)
    pattern = re.compile('\[([a-zA-Z0-9_\-]+)\]')
    # with gzip.open(opts.read_1, "rt") as f:
    with file_open(fastq) as f:
        for line in f:
            barcodes = pattern.findall(line.decode('utf-8'))

            if total < 5:
                assert len(barcodes) > 1, 'Barcode not found. BarcodeIndentification may not have been run on this fastq file!'
            #allow only a single missing barcode per read
            if barcodes.count('NOT_FOUND') < 2:
                #create a string to group barcodes by
                simp_barcode = '.'.join([convert_dict[i.split('_')[-1].strip('AB')] for i in barcodes])

                bc_container[simp_barcode].append(line.decode('utf-8').rstrip())
            else:
                skipped += 1
            next(f)
            next(f)
            next(f)
            total += 1

    print('Total reads:', total)
    print('Skipped reads (too many missing barcodes):', skipped)

    return bc_container


def parse_barcodes(fastq):
    '''Parse fastq barcodes into a dict
    '''
    skipped = 0
    total = 0
    bc_container = defaultdict(set)
    pattern = re.compile('\[([a-zA-Z0-9_\-]+)\]')

    with file_open(fastq) as f:
        for line in f:
            barcodes = pattern.findall(line.decode('utf-8'))
            # print(barcodes)
            if total < 5:
                assert len(barcodes) > 1, 'Barcode not found. BarcodeIndentification may not have been run on this fastq file!'

            #allow only a single missing barcode per read
            if barcodes.count('NOT_FOUND') < 2:
                #create a string to group barcodes by
                barcode = '|'.join(barcodes)
                bc_container[barcode].add(line.decode('utf-8').rstrip().split('::')[0])
            else:
                skipped += 1
            next(f)
            next(f)
            next(f)
            total += 1

    print('Total reads:', total)
    print('Skipped reads (too many missing barcodes):', skipped)

    return bc_container
#
# fastq = '/mnt/data/20190704_10_1000_spritezero/assembled/complex_10_S2_L001_R1_001_val_1.assembled_rv_bID.fq.gz'
# barcodes_10 = parse_barcodes(fastq)
# barcodes_10['R8-1_TermStag_bot_7|R7-2_bot_A1|R6-1_bot_B9|R5-2_bot_A3|R4-2_bot_A7|R3-1_bot_A5|R2-2_bot_B4|R1-1_bot_A2']
# bc_list = ham_adj_list(list(bc_container.keys())[0:1000])
# len(bc_container.keys())

def parse_clusters(path):
    '''Parse clusters file

    Args:
        path (str): Location of cluster file
    '''
    # path = '/mnt/data/20190704_10_1000_spritezero/assembled/10_unique.clusters'

    cluster_file = defaultdict(list)
    with file_open(path) as c_file:
        for line in c_file:
            barcode, *clusters = line.decode('utf-8').rstrip().split('\t')
            #Allow only one NOT_FOUND in barcodes
            if barcode.count('NOT_FOUND') < 2:
                cluster_file[barcode].extend(clusters)
    return cluster_file


def make_N_minus_1_barcodes(barcodes):
    '''For each barcode, make all N-1 possible barcodes
    '''
    # bc = 'R8-1_TermStag_bot_7|R7-2_bot_A1|R6-1_bot_B9|R5-2_bot_A3|R4-2_bot_A7|R3-1_bot_A5|R2-2_bot_B4|R1-1_bot_A2'
    bc_minus_n = defaultdict(set)
    for bc in barcodes.keys():
        bc_lst = bc.split('|')
        for i in range(len(bc_lst)):
            bc_minus_n[bc].add('|'.join(bc_lst[:i] + bc_lst[i+1:]))
    return bc_minus_n

def make_N_minus_X_barcodes(barcodes, X=1):
    '''For each barcode, make all N-X possible barcodes
    '''
    # bc = 'R8-1_TermStag_bot_7|R7-2_bot_A1|R6-1_bot_B9|R5-2_bot_A3|R4-2_bot_A7|R3-1_bot_A5|R2-2_bot_B4|R1-1_bot_A2'
    bc_minus_x = defaultdict(set)
    for bc in barcodes.keys():
        bc_lst = bc.split('|')
        items_to_remove = list(combinations(range(0,8),X))
        for y in items_to_remove:
            keep_items = [list(range(0,len(bc_lst)))[i] for i in range(0,len(bc_lst)) if i not in y]
            bc_minus_x[bc].add('|'.join([bc_lst[i] for i in keep_items]))
    return bc_minus_x


def get_barcodes2merge(barcodes, x=1):
    '''Get combination of barcodes that should be merged (i.e. )
    '''

    all_n1_bc = make_N_minus_X_barcodes(barcodes, x)

    #make inverted dict, this will give all the barcodes with N-1 differences
    invrt_all_n1_bc = defaultdict(set)
    for key, values in all_n1_bc.items():
        for i in values:
            invrt_all_n1_bc[i].add(key)

    #Merge all keys into barcode groups
    merged_barcodes = []
    for key, values in all_n1_bc.items():
        barcode_set = set()
        for item in values:
            for i in invrt_all_n1_bc.get(item):
                barcode_set.add(i)
        merged_barcodes.append(barcode_set)

    return merged_barcodes

def corrected_clusters(clusters, merged_barcodes, unique):
    '''Merge values of N-X clusters together
    '''
    if unique:
        clusters_out = defaultdict(set)
    else:
        clusters_out = defaultdict(list)

    for merge in merged_barcodes:
        new_bc = list(merge)[0] #just use the first for now
        for bc in merge:
            if unique:
                [clusters_out[new_bc].add(i) for i in clusters.get(bc)]
            else:
                clusters_out[new_bc].extend(clusters.get(bc))
    return clusters_out


def write_out(clusters, out_path):
    '''Write out corrected clusters file
    '''
    with open(out_path, 'w') as out_file:
        for key, values in clusters.items():
            out_file.write(key + '\t' + '\t'.join(values) + '\n')




def chunk_it(list_to_chunk, chunks):
    '''Split list into x smaller lists

    Args:
        list_to_chunk (list): A list of values to be chunked
        chunks (int): number of chunks to split the list into

    Return:
        A list of chunks (lists)
    '''
    avg = len(list_to_chunk) / float(chunks)
    out = []
    last = 0.0

    while last < len(list_to_chunk):
        out.append(list_to_chunk[int(last):int(last + avg)])
        last += avg

    return out


# @njit
def ham_adj_worker(umi_chunk, umis, threshold=1):
    '''Identify all umis within the hamming distance threshold (1 is best)

    Args:
        umi_chunk
        umis
        threshold (int): allowed hamming distance between barcodes
    '''

    adj_list = defaultdict(list)
    for i, umi in enumerate(umi_chunk):
        a1 = adj_list[umi]
        for j in range(i+1,len(umis)):
            umi2 = umis[j]
            if edit_distance(umi.encode('utf-8'), umi2.encode('utf-8')) <= threshold:
            # if distance.hamming(umi, umi2) <= threshold:
                adj_list[umi].append(umi2)

    return [adj_list]


# ham_adj_worker(list(bc_container.keys())[3], list(bc_container.keys())[0:1000])


def simplify_dict(in_dict):
    '''Simplify a dictionary of where all common keys and values are collapse into a single set

    Args:
        in_dict (dict): A dictionary where keys appear as values in other keys

    Example:
        Collapses a full dictionary
        'E.B.A.C.A.B.K.4': ['E.B.A.C.A.B.K.3','B.B.A.C.A.B.K.4','E.B.A.C.A.B.K.4']
        'E.B.A.C.A.B.K.3': ['E.B.A.C.A.B.K.3', 'E.B.A.C.A.B.K.4']
        into unique sets
        {'E.B.A.C.A.B.K.4','E.B.A.C.A.B.K.3','B.B.A.C.A.B.K.4','E.B.A.C.A.B.K.4'}

    Note:
        Messes up in_dict, need to make a deepcopy

    Return:
        List of sets
    '''

    k_v = [(set([k]), set(v)) for k, v in in_dict.items()]
    merged = 1
    while merged:
        merged = 0  #False
        results = []
        while k_v:
            common, rest = k_v[0], k_v[1:]
            k_v = []
            for x in rest:
                if x[1].isdisjoint(common[1]): #sets are disjoint only if their intersection is the empty set
                    k_v.append((x[0], x[1])) #don't share values with other sets
                else:
                    merged = 1 #True
                    common[1].update(x[1])
                    common[0].update(x[0])
            results.append(common)
        k_v = results
    #return only keys (shared keys)
    return [tp[0] for tp in k_v]



def get_corrected_qname(all_results_simp, bc_container):
    '''New qname for grouped barcodes

    Args:
        all_results_simp (list of sets): output from simplify_dict
        bc_container (dict): reads grouped by barcode
    Note:
        Based on grouped barcodes, change all reads to the same barcode.
        The barcode with the highest number of reads is used for all other reads.

    Example:
        Four barcodes in a group: ['C.G.K.B.G.C.D.4', 'G.G.K.B.G.C.D.4', 'G.G.K.G.G.C.D.4', 'B.G.K.B.G.C.D.4']
        'G.G.K.B.G.C.D.4' has the highest number of reads
        Therefore 'C.G.K.B.G.C.D.4', G.G.K.G.G.C.D.4', 'B.G.K.B.G.C.D.4' are changed 'G.G.K.B.G.C.D.4'

    Return:
        Dictionary of original read name (key), and the new read name (value)
    '''

    qname_correction = defaultdict()
    pattern = re.compile('\[([a-zA-Z0-9_\-]+)\]')
    #write out new barcodes
    #create a qname correction dictionary
    for barcodes in all_results_simp:
        #correct barcodes, assign the most frequent barcode to all other reads
        if len(barcodes) > 1:
            b_list = list(barcodes)
            read_counts = [len(bc_container[b]) for b in b_list]
            #if multiple barcodes have the same amount of reads, the barcode from the first one will be used
            barcode2use = read_counts.index(max(read_counts))
            out_barcode = pattern.findall(bc_container[b_list[barcode2use]][0])
            out_format = '[{0}]'.format(']['.join(out_barcode))
            for b in b_list:
                for read in bc_container[b]:
                    qname_correction[read] = read.split('::')[0] + '::' + out_format
        #if only a single barcode, just write out that barcode (no correction)
        else:
            for read in bc_container[list(barcodes)[0]]:
                qname_correction[read] = read

    return qname_correction

#
# read_counts = [len(bc_container[b]) for b in ['C.G.K.B.G.C.D.4', 'G.G.K.B.G.C.D.4', 'G.G.K.G.G.C.D.4', 'B.G.K.B.G.C.D.4']]
# barcode2use = read_counts.index(max(read_counts))
#
# barcodes = [bc_container[b] for b in ['C.G.K.B.G.C.D.4', 'G.G.K.B.G.C.D.4', 'G.G.K.G.G.C.D.4', 'B.G.K.B.G.C.D.4']]
#
# bc_container[list({'C.A.E.D.B.I.K.4'})[0]]



def write_corrected(qname_corrected, fq_path, out_path):
    '''Write fastq with new corrected barcodes

    Args:
        qname_corrected (dict): Dictionary with old (key) and new (value) qnames
        fq_path (str): input fastq path
        out_path (str): output fastq path
    '''

    if not out_path.endswith('.gz'):
        out_path = out_path + '.gz'

    with file_open(fq_path) as in_fq, \
    gzip.open(out_path, 'wt') as out_fq:
        for qname, seq, thrd, qual in tqdm(fastq_parse(in_fq)):
            try:
                new_qname = qname_corrected[qname]
                out_fq.write('\n'.join([new_qname, seq, thrd, qual]) + '\n')
            except KeyError:
                continue #skip reads with more than 1 missing barcode
                # raise Exception('Read missing')


def parse_args():

    parser = argparse.ArgumentParser(description='Split fastq based on dpm and rpm sequence')
    parser.add_argument('-c', '--clusters', dest='clusters', type=str, required=True,
                        help='Clusters file')
    parser.add_argument('-o', '--out', dest='out_path', type=str, required=True,
                        help='Output corrected clusters file')
    parser.add_argument('-x', '--xmimus', dest='xminus', type=int, default=1,
                        help='Number of errors allowed in barcode')

    # parser = argparse.ArgumentParser(description='Split fastq based on dpm and rpm sequence')
    # parser.add_argument('--r1', dest='read_1', type=str, required=True,
    #                     help='Fastq read 1')
    # parser.add_argument('--out', dest='out_path', type=str, required=True,
    #                     help='Output fastq.gz')
    # parser.add_argument('--threads', dest='threads', type=int, required=False,
    #                     default=1, help='Number of threads to use')

    opts = parser.parse_args()

    return opts


def main():

    opts = parse_args()

    # fastqgzfile = '/mnt/data/SPRITEzero_complexes/trimmed/assembled/25000_S1_L001_R1_001_val_1_assembled_rv_bID_10k.fastq.gz'
    # outfastqgzfile = '/mnt/data/SPIRITEzero_complexes/trimmed/assembled/25000_S1_L001_R1_001_val_1_assembled_rv_bID_corrected.fastq.gz'



    clusters = parse_clusters(opts.clusters)
    merge_barcodes = get_barcodes2merge(clusters, opts.xminus)

    out_clusters = corrected_clusters(clusters, merge_barcodes, unique=True)
    write_out(out_clusters, opts.out_path)

    #
    # clusters = parse_clusters('/mnt/data/20190704_10_1000_spritezero/assembled/10_unique.clusters')
    # merge_barcodes = get_barcodes2merge(clusters)
    #
    # out_clusters = corrected_clusters(clusters, merge_barcodes, unique=True)
    # write_out(out_clusters, '/mnt/data/20190704_10_1000_spritezero/assembled/10_unique_corrected.clusters')
    #
    #
    # clusters = parse_clusters('/mnt/data/20190704_10_1000_spritezero/assembled/1000_unique.clusters')
    # merge_barcodes = get_barcodes2merge(clusters)
    #
    # out_clusters = corrected_clusters(clusters, merge_barcodes, unique=True)
    # write_out(out_clusters, '/mnt/data/20190704_10_1000_spritezero/assembled/1000_unique_corrected.clusters')

    #
    #
    # bc_container = get_barcodes(opts.read_1)
    #
    # results = Parallel(n_jobs=opts.threads)(delayed(ham_adj_worker)(umi_chunk, list(bc_container.keys()), 1)
    #                               for umi_chunk in chunk_it(list(bc_container.keys()), opts.threads))
    #
    # results = ham_adj_worker(list(bc_container.keys(), list(bc_container.keys()), 1)
    #
    # #merge all dicts
    # all_results = {}
    # for d in results:
    #     all_results.update(d[0])
    #
    # assert len(all_results) == len(bc_container), 'Lenght post not same as pre'
    #
    # all_results_simp = simplify_dict(deepcopy(all_results))
    #
    # print('Clusters pre correction:', len(all_results))
    # print('Clusters post correction:', len(all_results_simp))
    #
    #
    # qname_corrected = get_corrected_qname(all_results_simp, bc_container)
    #
    # write_corrected(qname_corrected, opts.read_1, opts.out_path)



if __name__ == "__main__":
    main()
