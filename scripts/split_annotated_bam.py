
#%%
import pysam
import argparse
#%%
def parse_args():

    parser = argparse.ArgumentParser(description='Add tags to bam file')
    parser.add_argument('-i', '--input', dest='input_bam', type=str, required=True,
                        help='Input bam to split by featureCounts annotation')
    parser.add_argument('-ou', '--unanno', dest='output_unanno', type=str, required=True,
                        help='Unannotated bam out')
    parser.add_argument('-o', '--out_bam', dest='output_bam', type=str, required=True,
                        help='Annotated bam out')

    args = parser.parse_args()
    return args


# bam = "/mnt/data/RNA_DNA_SPRITE/featureCounts/PB49.RNAr.hisat2.mapq20.bam.featureCounts.bam"
# output_bam = "/mnt/data/RNA_DNA_SPRITE/featureCounts/PB49.RNAr.only.hisat2.mapq20.bam.featureCounts.bam"
# unanno_bam = "/mnt/data/RNA_DNA_SPRITE/featureCounts/PB49.RNAr.unanno.hisat2.mapq20.bam.featureCounts.bam"

def main():
    opts = parse_args()

    split_bams(opts.input_bam, opts.output_bam, opts.output_unanno)

#%%
def split_bams(bam, output_bam, unanno_bam):
    '''Split bams based on their featureCounts annotation presence

    Args:
        bam(str): input bam with all reads (annotated and unannotated)
        output_bam(str): reads only with featureCounts annotation
        unanno_bam(str): reads without annotation
    '''
    anno_reads = 0
    unanno_reads = 0
    with pysam.AlignmentFile(bam, 'rb') as in_bam, \
    pysam.AlignmentFile(output_bam, "wb", template=in_bam) as out_bam, \
    pysam.AlignmentFile(unanno_bam, "wb", template=in_bam) as un_bam:
        for read in in_bam.fetch(until_eof = True):
            #get featureCounts annotation
            if read.has_tag('XT'):
                out_bam.write(read)
                anno_reads += 1
            #only XS flag present if no feature identified
            elif read.has_tag('XS'):
                un_bam.write(read)
                unanno_reads += 1
            else:
                raise Exception('XS tag missing, was featureCounts run?')
                    
    print('Annotated reads:', anno_reads)
    print('Unannotated reads:', unanno_reads)
#%%
if __name__ == "__main__":
    main()