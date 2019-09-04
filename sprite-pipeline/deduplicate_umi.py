
#%%
from tqdm import tqdm
from collections import defaultdict, Counter
from babrahamlinkon.general import fastq_parse, file_open
from babrahamlinkon import deduplication_general
from joblib import Parallel, delayed
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import os


'''
Group sequences with the same UMI (9 nt UMI followed by CAAGTCA). Perform MSA to create a single
consensus sequence/quality and write out. Bundles (UMI groups) that have very different sequences
are filtered out.
'''

#%%
def make_bundle(fastq, umi_len, orientation='r1'):
    '''bundle reads using umi sequence
    '''

    reads_dict = defaultdict(dict)
    qual_dict = defaultdict(list)

    with file_open(fastq) as in_fq:
        for qname, seq, thrd, qual in fastq_parse(in_fq):

            #Get UMI from qname
            #9 nt UMI followed by CAAGTCA
            # s_qname, e_qname = qname.split(' ')

            if orientation == 'r1': #UMI on left
                umi = seq[:umi_len]
                trim_seq = seq[umi_len:]
                trim_qual = qual[umi_len:]
            elif orientation == 'r2': #UMI on right
                # seq = 'ACGTTGCGGGCATACCTTGGCGGGAGCGTTGTTGCGGATTGACGTAAGGCGGCGCTGATATTCCGCGCATCTGCTCGGCGTCTCAACAGCGCGCTTGTACGAGAGTCGGATGCTGACTTGTAGAGAATGTAAGGAA'
                umi = seq[-umi_len:]
                trim_seq = seq[:-umi_len]
                trim_qual = qual[:-umi_len]

            else:
                raise Exception('Orientation should be "r1" or "r2"')

            #create dictionary of sequence and quality...
            qual_dict[trim_seq].append(trim_qual)

            try:
                reads_dict[umi]['count'] += 1
                reads_dict[umi]['seq'].update([trim_seq]) #add all the seqs for consensus

            except KeyError:
                reads_dict[umi]['count'] = 1
                reads_dict[umi]['read'] = qname
                reads_dict[umi]['seq'] = Counter([trim_seq]) #add all the seqs for consensus

    #if same sequence has 2 quals take highest

    for k,v in qual_dict.items():
        if len(v) > 1:
            qual_lofls = [list(item) for item in v]
            qual_dict[k] = [deduplication_general.qual_highest(qual_lofls)]

    return (reads_dict, qual_dict)



#%%
class results():
    '''Holder for results
    '''
    def __init__(self):

                self.reads = []
                self.consensus_seqs = []
                self.final_umis = []
                self.umi_counts = []
                self.low_gt = 0
                self.corrected = 0
                self.low_gt_corrected = 0
                self.num_input = 0
                self.stats_pre_df_dict = {'UMI': [], 'counts': []}
                self.cons_diffs = defaultdict()
                self.cons_algn = defaultdict()
                self.cons_all = defaultdict(lambda: defaultdict())
                self.consensus_quals = []


def chunk_it(seq, num):
    '''
    Divide clusters into x chunks to process in parallel
    '''
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out




def reduce_clusters_worker(bundle, clusters, counts, qual_dict,
                           gt_threshold, cons_no_qual, mismtch, with_N):
    '''
    Parallel worker to reduce clusters into a single consensus sequence
    '''

    reads = []
    consensus_seqs = []
    consensus_quals = []
    final_umis = []
    umi_counts = []
    low_gt = 0
    corrected = 0
    low_gt_corrected = 0

    cons_diffs = defaultdict()
    cons_algn = defaultdict()


    for cluster in tqdm(clusters):
        umi_in_cluster = len(cluster)
        #Consensus for read loss filter
        out_dict = {umi:bundle[umi]['seq'] for umi in cluster}

        assert len(out_dict) > 0, 'No sequence from umi'

        if umi_in_cluster > 1: #contains umi with 1 error
            corrected += 1

        alignment, new_qual_dict = deduplication_general.kalign_msa(out_dict, qual_dict)
        gt_ratio, consensus_seq, consensus_qual, diffs_from_cons = \
            deduplication_general.read_loss(alignment, new_qual_dict, 
                                            differences=mismtch,
                                            no_msa=False, cons_no_qual=cons_no_qual,
                                            j_trim=0, with_N=with_N) #umi=umi_cons

        #keep record of the distance between sequence and consensus per UMI bases (clustered UMIs seperated by ,)
        if not isinstance(diffs_from_cons, int):
            cons_algn[','.join(cluster)] = ','.join(x for x in alignment)
            cons_diffs[','.join(cluster)] = ','.join(x for x in diffs_from_cons)
        else:
            cons_algn[','.join(cluster)] = diffs_from_cons
            cons_diffs[','.join(cluster)] = diffs_from_cons


        if gt_ratio >= gt_threshold:
            #Parent umi = highest count umi which account for the cluster
            parent_umi = deduplication_general.get_best_higher_counts(cluster, counts)
            reads.append(bundle[parent_umi]['read'])
            consensus_seqs.append(consensus_seq.replace('-', '')) #remove padding or indels from msa
            consensus_quals.append(consensus_qual.replace('#', '')) #should not have any # qual as preclean removed them

            final_umis.append(parent_umi)
            #Number of UMI's in the cluster (how many have been collapsed)
            umi_counts.append(sum([counts[x] for x in cluster]))
        else:
            low_gt += 1

            if umi_in_cluster > 1: #contains umi with 1 error and low ratio
                low_gt_corrected += 1

    return [reads, consensus_seqs, consensus_quals, final_umis, umi_counts, low_gt,
            corrected, low_gt_corrected, cons_diffs, cons_algn]





def deduplicate_bundle_parallel(bundle, qual_dict,
                    mismatch, stats, threads, gt_threshold,
                    cons_no_qual=False, with_N=False):
    '''
    Deduplicated reads in parralel using the reduce_clusters_worker
    '''

    dir_adj_results = []

    #do in parallel
    #reduce bundle
    umis = bundle.keys()

    bundle_results = results()
    #clusters need to be a list of sets
    clusters = []
    for key in bundle.keys():
        clusters.append(set([key]))

    print('Number of clusters to process:', len(clusters))

    list_of_clusters = chunk_it(clusters, threads)

    counts = {umi: bundle[umi]['count'] for umi in umis}


    #num_input
    bundle_results.num_input += sum([bundle[umi]['count'] for umi in bundle])

    # pre_average_distance = ''
    if stats:
        bundle_results.stats_pre_df_dict['UMI'].extend(bundle) #umi + read
        bundle_results.stats_pre_df_dict['counts'].extend([bundle[UMI]['count'] for UMI in bundle]) #umi counts
    # if stats:
    #     dir_adj_results[0][8]['UMI'].extend(bundle) #umi + read
    #     dir_adj_results[0][8]['counts'].extend([bundle[UMI]['count'] for UMI in bundle]) #umi counts

    dir_adj_results_lists = \
    Parallel(n_jobs=threads, backend='loky')(delayed(reduce_clusters_worker)(bundle, clusters, counts,
    qual_dict, gt_threshold, cons_no_qual, mismatch, with_N) for clusters in list_of_clusters)

    for reads_s, consensus_seqs_s, consensus_quals_s, final_umis_s, umi_counts_s, low_gt_s, \
    corrected_s, low_gt_corrected_s, cons_diffs_s, cons_algn_s in dir_adj_results_lists:

        bundle_results.reads.extend(reads_s)
        bundle_results.consensus_seqs.extend(consensus_seqs_s)
        bundle_results.consensus_quals.extend(consensus_quals_s)
        bundle_results.final_umis.extend(final_umis_s)
        bundle_results.umi_counts.extend(umi_counts_s)
        bundle_results.low_gt += low_gt_s
        bundle_results.corrected += corrected_s
        bundle_results.low_gt += low_gt_corrected_s
        for k,v in cons_diffs_s.items():
            try:
                bundle_results.cons_all[k]['diffs_from_cons'] += '_' + v

            except KeyError:
                bundle_results.cons_all[k]['diffs_from_cons'] = v


        for k,v in cons_algn_s.items():
            try:
                bundle_results.cons_all[k]['alignments'] += '_' + v
            except KeyError:
                bundle_results.cons_all[k]['alignments'] = v

    dir_adj_results.append(bundle_results)

    return dir_adj_results




def write_out_deduplicated(dir_adj_results, low_umi_out, out, stats, min_reads, pdf_out):
    '''
    Write out reads from result object produced from deduplication and
    return stats
    '''
    stats_pre_df_dict_all = {'UMI': [], 'counts': []}
    stats_post_df_dict = {'UMI': [], 'counts': []}
    pre_cluster_stats = []
    post_cluster_stats = []
    stats_cons_diffs = defaultdict(lambda: defaultdict())


    num_input_all, num_output = 0, 0

    low_gt_reads = 0
    corrected_reads = 0
    low_gt_corrected_reads = 0
    low_umi_count = 0


    print('Writing out')

    for bundle in range(len(dir_adj_results)):
        # reads, consensus, final_umis, umi_counts, low_gt, corrected, low_gt_corrected,\
        # topologies, nodes, num_input, stats_pre_df_dict, pre_average_distance = results
        # print('Unique:', len(set(dir_adj_results[bundle][0])), 'All:', len(dir_adj_results[bundle][0]))
        #umi_counts
        #in case no reads are output
        if dir_adj_results[bundle].umi_counts:
            #*-operator to unpack the arguments out of a list or tuple #dir_adj_results[bundle][3]
            labels, values = zip(*Counter(dir_adj_results[bundle].umi_counts).items())
            non_dedup_values = tuple(l*v for l, v in zip(labels, values))

            if min_reads != None:
                # cut_point = count_change(labels, non_dedup_values)

                cut_point = min_reads

                plt.figure()
                plt.bar(labels, non_dedup_values)

                plt.title(' Cut point: ' + str(cut_point), ha='center') #need to put name to know which bundle J is being processed
                my_plot = plt.axvline(cut_point, linestyle='dashed', linewidth=2).get_figure()
                pdf_out.savefig(my_plot)
                plt.close('all')

        num_input_all += dir_adj_results[bundle].num_input #dir_adj_results[bundle][7] #num_input

        #remove low umi counts 1-5
        indx = 0
        for count in dir_adj_results[bundle].umi_counts: #dir_adj_results[bundle][3]: #umi_counts
            #write out low quality UMI clusters
            if count <= cut_point:

                assert len(dir_adj_results[bundle].consensus_seqs[indx]) == len(dir_adj_results[bundle].consensus_quals[indx]), \
                'Consensus sequence not same length as consensus quality'
                # dir_adj_results[bundle][1]    dir_adj_results[bundle][10]

                low_umi_out.write(dir_adj_results[bundle].reads[indx].split(' ')[0] + '_' +
                                  dir_adj_results[bundle].final_umis[indx] + '_' + str(count) +'\n' +
                                  dir_adj_results[bundle].consensus_seqs[indx] + '\n' +
                                  '+' + '\n' +
                                  dir_adj_results[bundle].consensus_quals[indx] + '\n')
                low_umi_count += 1
            #write out good UMI clusters
            else:

                assert len(dir_adj_results[bundle].consensus_seqs[indx]) == len(dir_adj_results[bundle].consensus_quals[indx]), \
                'Consensus sequence not same length as consensus quality'
                out.write(dir_adj_results[bundle].reads[indx].split(' ')[0] + '_' +
                          dir_adj_results[bundle].final_umis[indx] + '_' + str(count) +'\n' +
                          dir_adj_results[bundle].consensus_seqs[indx] + '\n' +
                          '+' + '\n' +
                          dir_adj_results[bundle].consensus_quals[indx] + '\n')
                num_output += 1
            indx += 1


        if stats:
            # assert len(reads)  == len(consensus) == len(umi_counts), 'Reads, consensus and counts differ in length'

            ##################

            low_gt_reads += dir_adj_results[bundle].low_gt #dir_adj_results[bundle][4] #low_gt
            corrected_reads +=  dir_adj_results[bundle].corrected #dir_adj_results[bundle][5] #corrected
            low_gt_corrected_reads += dir_adj_results[bundle].low_gt_corrected #dir_adj_results[bundle][6] #low_gt_corrected

            # # collect pre-dudupe stats
            # stats_pre_df_dict['UMI'].extend(bundle) #umi + read
            # stats_pre_df_dict['counts'].extend([bundle[UMI]['count'] for UMI in bundle]) #umi counts
            #
            # pre_average_distance = get_average_umi_distance(bundle.keys()) #v_seq + umi
            stats_pre_df_dict_all['UMI'].extend(dir_adj_results[bundle].stats_pre_df_dict['UMI']) #stats_pre_df_dict umi + read #dir_adj_results[bundle][8]
            stats_pre_df_dict_all['counts'].extend(dir_adj_results[bundle].stats_pre_df_dict['counts']) #stats_pre_df_dict umi counts

            # pre_cluster_stats.append(pre_average_distance)

            #aggregate errors per cluster for each bundle
            cons_diffs = dir_adj_results[bundle].cons_all #dir_adj_results[bundle][9]
            for k,v in cons_diffs.items():
                try:
                    # print(stats_cons_diffs[k], '_', v)
                    stats_cons_diffs[k]['diffs_from_cons'] += '_' + v['diffs_from_cons']
                    stats_cons_diffs[k]['alignments'] += '_' + v['alignments']
                except KeyError:
                    stats_cons_diffs[k]['diffs_from_cons'] = v['diffs_from_cons']
                    stats_cons_diffs[k]['alignments'] = v['alignments']

            # collect post-dudupe stats
            #v_seq + umi
            # post_cluster_umis = [qname.split(' ')[-1] for qname in dir_adj_results[bundle][0]] #reads are just qnames
            stats_post_df_dict['UMI'].extend(dir_adj_results[bundle].final_umis) #final_umis #dir_adj_results[bundle][2]
            stats_post_df_dict['counts'].extend(dir_adj_results[bundle].umi_counts) #umi_counts #dir_adj_results[bundle][3]




    return [stats_pre_df_dict_all, stats_post_df_dict, pre_cluster_stats, post_cluster_stats,
    num_input_all, num_output, low_gt_reads, corrected_reads, low_gt_corrected_reads, low_umi_count,
    stats_cons_diffs]




def aggregate_stats_df(stats_df):
    ''' return a data frame with aggregated counts per UMI'''


    # agg_df_dict = {}

    total_counts = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=np.sum).T
    total_counts.rename(columns={'counts':'total_counts'}, inplace=True)
    median_counts = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=np.median).T
    median_counts.rename(columns={'counts':'median_counts'}, inplace=True)
    times_observed = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=len).T
    times_observed.rename(columns={'counts':'times_observed'}, inplace=True)

    agg_df = pd.concat([total_counts, median_counts, times_observed], axis=1)

    return agg_df


def cap_clusters(bundle):
    '''If UMI cluster has more than 1000 unique sequences, do not make consensus,
    but just discard

    Args:
        bundle(tuple): reads_dict and quals_dict
    '''
    umi_remove = []
    for k, v in bundle[0].items():
        if len(v['seq']) >= 1000:
            umi_remove.append(k)

    print('Will remove', len(umi_remove), 'clusters')

    for k in umi_remove:
        bundle[0].pop(k, None)
        bundle[1].pop(k, None)

    return(bundle)



def parse_args():

    parser = argparse.ArgumentParser(description='Deduplicate fastq (pre-barcodeID) using UMI in DPM/RPM')
    parser.add_argument('--r', dest='read', type=str, required=True,
                        help='Fastq file')
    parser.add_argument('--orient', dest='orient', type=str, default='r1',
                        help='Orientation of the Fastq; read_1 (use r1) or read_2 reverse complement (use r2)')
    parser.add_argument('--umi_len', dest='umi_len', type=int, default=16,
                        help='Length of UMI sequence')
    parser.add_argument('--threads', dest='threads', type=int, default=1,
                        help='Number of threads to use for parallel computation')
    parser.add_argument('--count_cut', dest='cut', type=int, default=0,
                        help='Number of reads per cluster in order to retain that cluster')
    opts = parser.parse_args()

    return opts





def main():

    #argparse
    opts = parse_args()

    # set up arrays to hold stats data
    stats_pre_df_dict_all = {'UMI': [], 'counts': []}
    stats_post_df_dict_all = {'UMI': [], 'counts': []}
    pre_cluster_stats_all = []
    post_cluster_stats_all = []


    num_input_an1, num_output_an1 = 0, 0
    num_input_an2, num_output_an2 = 0, 0
    # line_fmt = "@{0!s}\n{1!s}\n+\n{2!s}\n"
    low_gt_reads_an1, low_gt_reads_an2 = 0, 0
    corrected_reads_an1, corrected_reads_an2 = 0, 0
    low_gt_corrected_reads_an1, low_gt_corrected_reads_an2 = 0, 0
    low_umi_count_an1, low_umi_count_an2 = 0, 0


    # input_fq = '/home/chovanec/Documents/SPIRITEzero_complexes/trimmed/assembled/25000_S1_L001_R1_001_val_1_assembled.fastq'
    # out_name = '/home/chovanec/Documents/SPIRITEzero_complexes/trimmed/assembled/25000_S1_L001_R1_001_val_1_assembled_umi.fastq'
    # low_name = '/home/chovanec/Documents/SPIRITEzero_complexes/trimmed/assembled/25000_S1_L001_R1_001_val_1_assembled_low_count_umi.fastq'
    # out_prefix = '/home/chovanec/Documents/SPIRITEzero_complexes/trimmed/assembled/25000_S1_L001_R1_001_val_1_assembled_umi'

    input_fq = os.path.abspath(opts.read)
    print('Processing:', input_fq)
    out_name = input_fq.split('.')[0] + '_umi.fastq'
    low_name = input_fq.split('.')[0] + '_low_count_umi.fastq'
    out_prefix = input_fq.split('.')[0] + '_umi'

    read_bundles_all = make_bundle(input_fq, opts.umi_len, opts.orient)

    #examine bundles to check for extremely large umi clusters that are causing memory leaks
    # test_bundle = make_bundle('/mnt/data/20190621_PB_PC_SQ/assembled/PB-PC_S2_L001_R1_001_val_1.assembled.fastq.gz', umi_len=16, orientation='r2')
    largest_bundle = max([len(i['seq']) for i in read_bundles_all[0].values()])
    print('Largest cluster:', largest_bundle)


    #cap cluster sizes at 1000 to avoid memory leak warnings
    read_bundles = cap_clusters(read_bundles_all)

    deduplication_results = deduplicate_bundle_parallel(read_bundles[0], read_bundles[1],
                                                        mismatch=5, gt_threshold=1,
                                                        stats=True, threads=opts.threads,
                                                        cons_no_qual=False,
                                                        with_N=False)

    with open(out_name, 'w') as fq_out, open(low_name, 'w') as low_umi_out, \
    PdfPages(out_prefix + '_histogram.pdf') as pdf:

        stats_pre_df_dict, stats_post_df_dict, pre_cluster_stats, post_cluster_stats, \
        num_input, num_output, low_gt_reads, corrected_reads, \
        low_gt_corrected_reads, low_umi_count, stats_cons_diffs=\
        write_out_deduplicated(deduplication_results, low_umi_out, fq_out, stats=True, min_reads=opts.cut,
                               pdf_out=pdf)

    num_input_an1 += num_input
    num_output_an1 += num_output
    low_gt_reads_an1 += low_gt_reads
    corrected_reads_an1 += corrected_reads
    low_gt_corrected_reads_an1 += low_gt_corrected_reads
    low_umi_count_an1 += low_umi_count

    stats_pre_df_dict_all.update(stats_pre_df_dict)
    stats_post_df_dict_all.update(stats_post_df_dict)

    pre_cluster_stats_all.extend(pre_cluster_stats)
    post_cluster_stats_all.extend(post_cluster_stats)

    stats_pre_df = pd.DataFrame(stats_pre_df_dict_all)
    stats_post_df = pd.DataFrame(stats_post_df_dict_all)

    # print(pd.DataFrame.from_dict(Counter(stats_post_df_dict['counts']), orient='index').reset_index())

    # generate histograms of counts per UMI at each position
    UMI_counts_df_pre = pd.DataFrame(stats_pre_df.pivot_table(
        columns=stats_pre_df['counts'], values='counts', aggfunc=len)).T #not sure why it need to be transposed now when it didnt before!

    UMI_counts_df_post = pd.DataFrame(stats_post_df.pivot_table(
        columns=stats_post_df['counts'], values='counts', aggfunc=len)).T

    # UMI_counts_df_pre.columns = ['instances']
    # UMI_counts_df_post.columns = ['instances']
    UMI_counts_df_pre.rename(columns={'counts':'instances'}, inplace=True)
    UMI_counts_df_post.rename(columns={'counts':'instances'}, inplace=True)

    # print(stats_pre_df.pivot_table(columns="UMI", values="counts", aggfunc=np.sum))

    UMI_counts_df = pd.merge(UMI_counts_df_pre, UMI_counts_df_post,
                             how='outer', left_index=True, right_index=True,
                             sort=True, suffixes=['_pre', '_post'])

    UMI_counts_df = UMI_counts_df.fillna(0).astype(int)

    UMI_counts_df.to_csv(out_prefix + '_per_umi_per_position.tsv', sep='\t')

    ##########################

     # aggregate stats pre/post per UMI
    agg_pre_df = aggregate_stats_df(stats_pre_df)
    agg_post_df = aggregate_stats_df(stats_post_df)

    agg_df = pd.merge(agg_pre_df, agg_post_df, how='left',
                      left_index=True, right_index=True,
                      sort=True, suffixes=['_pre', '_post'])

    agg_df = agg_df.fillna(0).astype(int)

    # stats_consensus_difference = pd.DataFrame(stats_cons_diffs, index=[0])
    stats_consensus_difference = pd.DataFrame(stats_cons_diffs)
    stats_consensus_difference = stats_consensus_difference.T
    stats_consensus_difference.columns = ['Alignments', 'Consensus_differences']

    # stats_consensus_difference['UMI'] = stats_consensus_difference.index

    agg_df = pd.merge(agg_df, stats_consensus_difference, how='left',
                      left_index=True, right_index=True,
                      sort=True,)

    agg_df.to_csv(out_prefix + '_per_umi.tsv', sep="\t")


    print('Number of input reads:', num_input_an1)
    print('Number of output reads:', num_output_an1)
    print('Number of clusters with low ratio discarded:' + str(low_gt_reads_an1))
    print('Number of low UMI count groups:', low_umi_count_an1)



if __name__ == "__main__":
    main()


#%%
