from itertools import combinations
from cluster import parse_cluster
from collections import defaultdict
import pandas as pd
from natsort import natsorted, index_natsorted, order_by_index
import pyranges as pr
import argparse

# d = {'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 8, 5], 'End': [6, 9, 7], 
#      'Score': [0.1, 5, 3.14], 'Strand': ['+', '+', '-']}
# gr = pr.from_dict(d)
# gr.to_rle()

def parse_arguments():
    parser = argparse.ArgumentParser(description =
            "Clusters to other format conversion")
    parser.add_argument('-i', '--input', action = 'store', metavar = 'FILE',
                        help = 'Input cluster file')
    parser.add_argument('-o', '--output', action = 'store', metavar = 'FILE',
                        help = 'Output file')
    parser.add_argument('--max_cluster_size', metavar = 'MAX', type = int,
                        action = 'store', default = 1000,
                        help = "Skip read-clusters with more reads than MAX. (default 1000)")
    parser.add_argument('--min_cluster_size', metavar = 'MIN', type = int,
                        action = 'store', default = 2,
                        help = "Skip read-clusters with fewer reads than MIN. (default 2)") 
    parser.add_argument('--normalise', action = 'store_true',
                        help = 'Should normalisation of cluster size be applied')
    parser.add_argument('--format', metavar = 'FORMAT', type = str,
                        action = 'store', default = 'sfws',
                        choices=['sfws', 'bed', 'h5'],
                        help = "What format to output") 
                       
    return parser.parse_args()


def main():

    args = parse_arguments()

    clusters = parse_cluster(args.input)

    if args.format == 'sfws':
        convert_clusters(clusters, args.min_cluster_size, args.max_cluster_size, args.output, args.normalise)

    if args.format == 'bed':
        c_pyr = cluster2pyranges(clusters, args.min_cluster_size, args.max_cluster_size, args.normalise)
        c_pyr.to_bed(args.output)

    if args.format == 'h5':
        c_pyr = cluster2pyranges(clusters, args.min_cluster_size, args.max_cluster_size, args.normalise)
        c_pyr.to_bed(args.output)

        clusters_df = c_pyr.df
        dpm_df = clusters_df[clusters_df['read_type']=='DPM']
        dpm_df = dpm_df.drop(['read_type', 'exon', 'intron', 'repeat'], axis=1)
        rpm_df = clusters_df[clusters_df['read_type']=='RPM']
        rpm_df = rpm_df.drop(['read_type', 'Score'], axis=1)
        store = pd.HDFStore(args.output, 'w', complevel=9, complib='lzo')
        # store = pd.HDFStore('/mnt/data/RNA_DNA_SPRITE/20200315_RNA_DNA_Mario/8226.comboall.h5', 'w',
        #                     complevel=9, complib='lzo')
        store.append('DPM', dpm_df, format='table', append=True, data_columns=True)
        store.append('RPM', rpm_df, format='table', append=True, data_columns=True)

        store.close()

    # clusters = parse_cluster('/mnt/data/RNA_DNA_SPRITE/20200315_RNA_DNA_Mario/8226.comboall.clusters.gz')
    # convert_clusters(clusters, 2,1000, '/mnt/data/RNA_DNA_SPRITE/20200315_RNA_DNA_Mario/8226.comboall.DPM.sfws.txt',
    #                  normalise=False)

    # c_pyr = cluster2pyranges(clusters, 2, 1000)
    # c_pyr.to_bed('/mnt/data/RNA_DNA_SPRITE/20200315_RNA_DNA_Mario/H1.comboall.bed')




def convert_clusters(clusters, cluster_size_min, cluster_size_max, out_path, normalise):
    '''
    Args:
        clusters(Object): Clusters object
        cluster_size_min(int): minimum size of cluster (more than or equal to int)
        cluster_size_max(int): maximum size of cluster (less than or equal to int)
        out_path(str): Out file path
    '''
    assert cluster_size_min > 1, 'Minimum cluster size needs to be > 1'

    all_clusters = []
    for barcode, cluster in clusters.get_items():
        cs = cluster.size('DPM')
        if cs >= cluster_size_min and cs <= cluster_size_max:
            all_clusters.extend(cluster2sfws(cluster, 'DPM', normalise))

    column_names=['str1', 'chr1', 'pos1', 'frag1', 'str2', 'chr2', 'pos2','frag2', 'score']
    df = pd.DataFrame(all_clusters, columns=column_names) 
    df_out = df.reindex(index=order_by_index(df.index, index_natsorted(zip(df.chr1, df.chr2, df.pos1, df.pos2))))
    df_out.to_csv(out_path, sep=' ', index=False, header=False)



def cluster2sfws(cluster, read_type, normalise=True):
    '''Convert a cluster class object (a single cluster)
    to a dictionary in the sfws format (Juicer tools Pre format with score)
    
    Note:
        Juicer short format with score (sfws)
        A whitespace separated file that contains, on each line
            <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <score>
        https://github.com/aidenlab/juicer/wiki/Pre

        IMPORTANT NOTE pre throws away reads that map to the same restriction fragment. 
        If you use dummy numbers for the frag field, be sure they are different for the 
        different read ends; that is, <frag1> should be 0 and <frag2> should be 1.

        str = strand (0 for forward, anything else for reverse)

    Args:
        cluster(Cluster): A single Cluster object that holds all the position (reads)
        read_type(str): RPM or DPM
    '''

    cluster_pos = []
    for position in cluster:
        if position._type == read_type:
            cluster_pos.append(('0' if position._strand == '+' else '1',
                                position._chromosome,
                                position._start_coordinate))
    
    pairs = list(combinations(cluster_pos, 2))
    if normalise:
        score = 2.0 / len(cluster_pos)
    else:
        score = 1
    out = []
    for a, b in pairs:
        # chr1 > chr2 order
        a, b = order_by_index([a, b], index_natsorted([a[1], b[1]]))
        out.append([*a, 0, *b, 1, score]) 
    return out



def cluster2pyranges(clusters, cluster_size_min, cluster_size_max, normalise=True):
    '''Convert cluster format to pyranges (BED)

    Notes:
        chrom, start, end, strand, cluster, exon, intron, repeat, read_type

    Args:
        cluster(Cluster): A single Cluster object that holds all the position (reads)
    '''
    chrom = list()
    start = list()
    end = list()
    strand = list()
    scores = list()
    barcodes = list()
    exon = list()
    intron = list()
    repeat = list()
    read_type = list()

    unq_barcodes = set()
    for br, cluster in clusters.get_items():
        #make sure barcodes are unique
        count = 0
        while br in unq_barcodes:
            br = br + '_' + str(count)
            count += 1
        else:
            unq_barcodes.add(br)

        cs = cluster.size()
        csd = cluster.size('DPM')
        if csd == 0:
            print(cs, csd)
        if normalise:
            score = 2/csd
        else:
            score = 1

        if cs >= cluster_size_min and cs <= cluster_size_max:
            for position in cluster:
                chromosome = position._chromosome if position._chromosome.startswith('chr') else 'custom'
                if chrom == 'custom':
                    assert position._type=='RPM', 'DPM is not aligned to custom'

                chrom.append(chromosome)
                start.append(position._start_coordinate)
                end.append(position._end_coordinate)
                strand.append(position._strand)              
                barcodes.append(br)
                if position._type == 'RPM':
                    scores.append(1)
                elif position._type == 'DPM':
                    scores.append(score)
                read_type.append(position._type)
                #annotation
                features = position._feature.split(';')
                fs = classify_feature(features)
                exon.append(fs.get('exon', 'NA'))
                intron.append(fs.get('intron', 'NA'))
                repeat.append(fs.get('repeat', 'NA'))

    #Convert to pyranges
    r_df = pd.DataFrame({'Chromosome': chrom, 'Start': start, 'End': end, 
                        'Strand': strand, 'Name': barcodes, 'Score': score, 'exon': exon, 
                        'intron': intron, 'repeat': repeat, 'read_type': read_type})

    gr = pr.PyRanges(r_df)
    
    return gr


def classify_feature(features):
    '''classify features

    Args:
        features(list): List of features to classify
    '''

    c_feature = {'exon':'NA', 'intron':'NA', 'repeat':'NA'}

    for f in features:
        if f.endswith('exon'):
            c_feature['exon'] = f
        elif f.endswith('intron'):
            c_feature['intron'] = f
        elif f.endswith('none'):
            continue
        elif not f.endswith('exon|intron|none') or f.endswith('repeat'):
            if len(f) > 0:
                c_feature['repeat'] = f

    return c_feature




if __name__ == "__main__":
    main()
    