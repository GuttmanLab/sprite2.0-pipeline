import argparse
import os
import gzip

def parse_args():

    parser = argparse.ArgumentParser(description='Split fastq based on cell barcode')
    parser.add_argument('-r1', '--read_1', dest='r1', type=str, required=True,
                        help='Read 1')
    parser.add_argument('-r2', '--read_2', dest='r2', type=str, required=True,
                        help='Read 2')
    parser.add_argument('-psm44', '--psm44_dir', dest='psm44_dir', type=str, required=True,
                        help='Output directory for psm44')
    parser.add_argument('-bsps', '--bsps_dir', dest='bsps_dir', type=str, required=True,
                        help='Output directory for bsps')
    args = parser.parse_args()
    return args


def main():

    opts = parse_args()

    psm44_r1_out = os.path.join(opts.psm44_dir + os.path.basename(opts.r1))
    psm44_r2_out = os.path.join(opts.psm44_dir + os.path.basename(opts.r2))
    bsps_r1_out = os.path.join(opts.bsps_dir + os.path.basename(opts.r1))
    bsps_r2_out = os.path.join(opts.bsps_dir + os.path.basename(opts.r2))

    filter_paired_end(opts.r1, opts.r2, ['PSM44', 'BSPS'], psm44_r1_out, psm44_r2_out,
                      bsps_r1_out, bsps_r2_out)



#%%
def filter_paired_end(fastq1, fastq2, split_list, ofq1_val1, ofq2_val2, ofq3_val1, ofq4_val2):
    '''Loop through both fastq's and separate cell specific records

    Args:
        fastq1(str): Path to fastq1
        fastq2(str): Path to fastq2
        split_list(list): A list of values on which to split the fastq into two files
        ofq1_val1(str): Output path of records in split list R1
        ofq2_val2(str): Output path of records in split list R2
        ofq3_val1(str): Output path of records not in split list R1
        ofq4_val2(str): Output path of records not in split list R2
    '''
    assert len(split_list) == 2, 'Number of split values needs to be 2'

    reads_dropped = 0
    with file_open(fastq1) as f1, file_open(fastq2) as f2, \
        gzip.open(ofq1_val1, 'wt') as psm1, \
        gzip.open(ofq2_val2, 'wt') as psm2, \
        gzip.open(ofq3_val1, 'wt') as bsps1, \
        gzip.open(ofq4_val2, 'wt') as bsps2:
        for fq1_line, fq2_line in zip(fastq_parse(f1), fastq_parse(f2)):
            qname1, seq1, thrd1, qual1 = fq1_line
            if split_list[0] in qname1:
                psm1.write('\n'.join(fq1_line) + '\n')
                psm2.write('\n'.join(fq2_line) + '\n')
            elif split_list[1] in qname1:
                bsps1.write('\n'.join(fq1_line) + '\n')
                bsps2.write('\n'.join(fq2_line) + '\n')
            else:
                #drop read
                reads_dropped += 1

    print(reads_dropped)

#%%
# r1 = '/mnt/data/20190711_PB_PC_sz_quench/no_quench_S1_L001_R1_001_val_1.fq.gz'
# r2 = '/mnt/data/20190711_PB_PC_sz_quench/no_quench_S1_L001_R2_001_val_2.fq.gz'

# count = 0
# with file_open(r1) as f1, file_open(r2) as f2:
#     for line1, line2 in zip(fastq_parse(f1), fastq_parse(f2)):
#         qname1, seq1, thrd1, qual1 = line1
#         if count < 100:
#             print(qname1)
#             count += 1
#         else:
#             break


#%%
def file_open(filename):
    """
    Open as normal or as gzip
    Faster using zcat?
    """
    #does file exist?
    f = open(filename,'rb')
    if (f.read(2) == b'\x1f\x8b'): #compressed alsways start with these two bytes
        f.seek(0) #return to start of file
        return gzip.GzipFile(fileobj=f, mode='rb')
    else:
        f.seek(0)
        return f


#%%
def fastq_parse(fp):
    """
    Parse fastq file.
    """
    linecount = 0
    name, seq, thrd, qual = [None] * 4
    for line in fp:

        linecount += 1
        if linecount % 4 == 1:
            try:
                name = line.decode('UTF-8').rstrip()
            except AttributeError:
                name = line.rstrip()
            assert name.startswith('@'),\
                   "ERROR: The 1st line in fastq element does not start with '@'.\n\
                   Please check FastQ file near line number %s" % (linecount)
        elif linecount % 4 == 2:
            try:
                seq = line.decode('UTF-8').rstrip()
            except AttributeError:
                seq = line.rstrip()
        elif linecount % 4 == 3:
            try:
                thrd = line.decode('UTF-8').rstrip()
            except AttributeError:
                thrd = line.rstrip()
            assert thrd.startswith('+'),\
                   "ERROR: The 3st line in fastq element does not start with '+'.\n\
                   Please check FastQ file near line number %s" % (linecount)
        elif linecount % 4 == 0:
            try:
                qual = line.decode('UTF-8').rstrip()
            except AttributeError:
                qual = line.rstrip()
            assert len(seq) == len(qual),\
                    "ERROR: The length of Sequence and Quality aren't equal.\n\
                    Please check FastQ file near line number %s" % (linecount)

            yield name, seq, thrd, qual,
            name, seq, thrd, qual = [None] * 4


if __name__ == "__main__":
    main()

#%%
