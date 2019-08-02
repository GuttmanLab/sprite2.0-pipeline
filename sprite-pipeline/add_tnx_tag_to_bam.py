
import pysam
import argparse

# my_bam = '/home/chovanec/2p5-2.RNAr.Aligned.sortedByCoord.out.bam'
# output_bam = '/home/chovanec/2p5-2.RNAr.Aligned.sortedByCoord.out.mod.bam'



def parse_args():

    parser = argparse.ArgumentParser(description='Add tags to bam file')
    parser.add_argument('-i', '--input_bam', dest='input_bam', type=str, required=True,
                        help='BAM path to which to add XS tag from chrom field')
    parser.add_argument('-o', '--out_bam', dest='output_bam', type=str, required=True,
                        help='Out BAM path')

    args = parser.parse_args()
    return args

def main():
    opts = parse_args()

    add_tag(opts.input_bam, opts.output_bam)


def add_tag(input_bam, output_bam):
    count = 0
    skipped = 0
    with pysam.AlignmentFile(input_bam, "rb") as in_bam, \
    pysam.AlignmentFile(output_bam, "wb", template=in_bam) as out_bam:

        for read in in_bam.fetch(until_eof = True):
            count += 1
            try:
                read.tags += [('XS', read.reference_name)]
                out_bam.write(read)
            except KeyError:
                skipped += 1

    print('Total reads:', count)
    print('Reads with an error not written out:', skipped)

if __name__ == "__main__":
    main()
