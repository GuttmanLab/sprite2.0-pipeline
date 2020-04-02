import argparse
from collections import defaultdict


def main():

    args = parse_arguments()

    cluster_dict = combine_clusters(args.input, args.ignore)
    write_cluster(cluster_dict, args.output)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Merge multiple cluster files')
    parser.add_argument('-i', '--input',
                        metavar = "FILE",
                        action = "store",
                        required=True,
                        nargs='+',
                        help = "The input cluster file(s).")
    parser.add_argument('-o', '--output',
                        metavar = "FILE",
                        action = "store",
                        required=True,
                        help = "The output merged cluster file.")
    parser.add_argument('--ignore',
                        action="store_true",
                        required=False,
                        default=False, 
                        help = "When merging, will ignore file name in barcode")
    
    return parser.parse_args()


# cluster_paths = ['/mnt/data/RNA_DNA_SPRITE/HEKRNADNASPRITE/workup/clusters/test_1.cluster', 
                #  '/mnt/data/RNA_DNA_SPRITE/HEKRNADNASPRITE/workup/clusters/test_2.cluster']

def combine_clusters(cluster_paths, ignore_file):
    '''Combine multiple cluster files into a single file

    Args:
        cluster_paths(list): List of cluster file paths
        ignore_file(logical): should name in barcode be ignored
    '''
    assert len(cluster_paths) > 1, 'Only a single cluster file provided'

    in_clusters = 0

    cluster_merge = defaultdict(str)

    for cluster_path in cluster_paths:
        with open(cluster_path) as in_cluster:
            for line in in_cluster:
                in_clusters += 1
                barcode, *reads = line.rstrip().split('\t')
                if ignore_file:
                    *barcode, file_name = barcode.split('.')
                    barcode = '.'.join(barcode)
                cluster_merge[barcode] += '\t' + '\t'.join(reads)

    print('In clusters:', in_clusters)
    print('Clusters after merging:', len(list(cluster_merge.keys())))

    return cluster_merge



def write_cluster(cluster_dict, out_path):
    '''From a dictionary of barcodes (keys) and reads (values) write out file
    
    Args:
        cluster_dict(dict): A dictionary with barcodes (keys) and reads (values)
        out_path(str): Output path
    '''

    clusters_written = 0    
    with open(out_path, 'w') as out:
        for barcode, reads in cluster_dict.items():
            clusters_written += 1
            out.write(barcode + reads + '\n')

    print('Clusters written out:', clusters_written)


if __name__ == "__main__":
    main()


# barcode_to_positions = defaultdict(list)

# with open(sys.argv[1], 'r') as f:
#     for line in f:
#         fields = line.rstrip().split()
#         barcode_to_positions[fields[0]] = fields[1:]

# with open(sys.argv[2], 'r') as f:
#     for line in f:
#         fields = line.rstrip().split()
#         print fields[0] + "\t" + "\t".join(set(barcode_to_positions[fields[0]] + fields[1:]))
#         del barcode_to_positions[fields[0]]

# for barcode, positions in barcode_to_positions.iteritems():
#     print barcode + "\t" + "\t".join(positions)


