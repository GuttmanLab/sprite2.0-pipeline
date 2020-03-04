import cluster as c
from collections import Counter
import argparse
import os


def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Create a consensus fasta file from dfam')
    parser.add_argument('-i', '--input',
                        type = str,
                        action = "store",
                        required=True,
                        nargs="+",
                        help = "Input files")
    parser.add_argument('-o', '--output',
                        type = str,
                        action = 'store',
                        required=True,
                        help = "Output file")

    return parser.parse_args()


def main():

    args = parse_arguments()

    with open(args.output, 'w') as out:
        #write header
        out.write('\t'.join(['File_name', 'Total_DPM', 'Total_RPM', 
                             'DPM_only', 'RPM_only']) + '\n')
        for f in args.input:

            file_name = os.path.basename(f)

            clusters = c.parse_cluster(f)

            total_counts = Counter()
            RPM_only = Counter()
            DPM_only = Counter()
            for cluster in clusters:
                cl_count = cluster.count_type()
                if cl_count['DPM'] == 0:
                    RPM_only += cl_count
                elif cl_count['RPM'] == 0:
                    DPM_only += cl_count
                total_counts += cl_count
            
            out.write('\t'.join([file_name, str(total_counts['DPM']), 
                                 str(total_counts['RPM']),
                                 str(DPM_only['DPM']),
                                 str(RPM_only['RPM'])]) + '\n')
        

if __name__ == "__main__":
    main()