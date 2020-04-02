#%%
import pyBigWig
import sys
import assembly
import pyranges as pr
import tqdm
import argparse
import re 
import pandas as pd
# import pyranges_db as pr_db



def parse_arguments():
    parser = argparse.ArgumentParser(description =
            "Bedgraph to bigwig conversion")
    parser.add_argument('-i', '--input', action = 'store', metavar = 'FILE',
                        help = 'Input bedgraph file')
    parser.add_argument('-o', '--output', action = 'store', metavar = 'FILE',
                        help = 'Output bigwig file')
    parser.add_argument('--assembly', metavar = "ASSEMBLY",
                        action = 'store', dest='assembly', default='mm10',
                        choices = ["mm9", "mm10", "hg19", "hg38"],
                        help = 'The genome assembly. (default mm10)')
    return parser.parse_args()

#%%     

def main():

    args = parse_arguments()

    chrom_sizes = assembly.build(args.assembly, 1)._chromsizes

    #chromsizes to pyranges
    chroms = []
    start = []
    end = []
    for k, v in chrom_sizes.items():
        chroms.append(k)
        start.append(0)
        end.append(v)

    chromsize_gr = pr.PyRanges(chromosomes=chroms, starts=start, ends=end)

    # bedgraph_path = '/mnt/data/RNA_DNA_SPRITE/U1b2.100bps.bedgraph'
    bedgraph_path = args.input

    bg = read_bedgraph(bedgraph_path, chrom_sizes)

    write_bigwig(bg, args.output, chrom_sizes, chromsize_gr)


#%%

class bedgraph:
    '''Bedgraph holder
    '''
    def __init__(self, chrom, start, end, score):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.score = score
        self.length = len(self.chrom)

    def convert2pyranges(self):
        r_df = pd.DataFrame({'Chromosome': self.chrom, 'Start': self.start, 
                      'End': self.end, 'Score': self.score})
        gr = pr.PyRanges(r_df)
        return gr

#%%
def read_bedgraph(bg_path, chrom_sizes):
    '''Read bedgraph file and correct end of chromosome coordinates
    
    Args:
        bg_path (str): File path of bedgraph 
        chrom_sizes (Ordered_dict): chromosome sizes from assembly.py
    
    '''

    chrom = list()
    start_coord = list()
    end_coord = list()
    value = list()
    with open(bg_path, 'r') as in_bg:
        for line in in_bg:
            chrm, strt, ed, val = line.rstrip('\n').split('\t')
            
            chrom.append(chrm)
            start_coord.append(int(strt))
            end_coord.append(int(ed))
            value.append(float(val))

    bg_out = bedgraph(chrom, start_coord, end_coord, value)
    return bg_out


#%%
def sorted_nicely(l): 
    ''' Sort the given iterable in the way that humans expect.
    https://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python
    '''
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

#%%
# chromsizes = pr_db.ucsc.chromosome_sizes("mm10")
# gr.to_bigwig("chipseq.bw", chromsizes)

def write_bigwig(bg_class, out_path, chrom_sizes, chromsize_gr, one_based=True):
    '''Write out bigwig
    
    '''
    gr = bg_class.convert2pyranges()
    # gr_reduced = gr.merge(count=True, slack=-1)
    #clip coordinates to the end of chromosomes
    gr_clipped = pr.gf.genome_bounds(gr, chromsize_gr, clip=True)

    subset = ['Chromosome', 'Start', 'End', 'Score']
    gr_clipped = gr_clipped[subset].unstrand()
    df = gr_clipped.sort(strand=False).df
    # df.to_csv('/mnt/data/RNA_DNA_SPRITE/U1_Malat1_bedgraphs/test.bg', sep='\t', index=False)

    # gr_clipped.to_bigwig(out_path, chromsize_gr)

    #header chromosome needs to be the same as order of records
    unique_chromosomes = list(set(df['Chromosome']))
    sorted_chroms = sorted_nicely(unique_chromosomes)
    header = [(c, int(chrom_sizes[c])) for c in sorted_chroms]

    bw = pyBigWig.open(out_path, "w")
    bw.addHeader(header)

    gr = gr[subset].unstrand()
    df = gr.sort(strand=False).df

    chrom = list(df['Chromosome'])
    start = list(df['Start'])
    end = list(df['End'])
    value = list(df['Score'])

    previous_end = 0
    previous_start = 0
    current_chrom = 'chr1'
    for i in range(len(chrom)):
        if current_chrom == chrom[i]:
            assert start[i] >= previous_start, f'Error not sorted start {previous_start} {chrom[i]} {start[i]}'
            assert end[i] >= previous_end, f'Error not sorted end {end[i]}'
            
        else:
            current_chrom = chrom[i]
        
        previous_start = start[i]
        previous_end = end[i]   


    bw.addEntries(chrom, start, ends=end, values=value, validate=True)

    bw.close()




if __name__ == "__main__":
    main()
    



#%%
# def resolve_conflicts(bg_class, chromsize_gr):
#     '''
#     Make sure no overlapping coordinate are present in bedgraph
#     Make sure coordinates don't run past chromosome end

#     Args:
#         bg_class(object): bg class object
#         chromsize_gr(pyranges object): pyranges with chromosome sizes
#     '''
#     gr = bg_class.convert2pyranges()
#     gr_reduced = gr.merge(count=True, slack=-1)
#     #clip coordinates to the end of chromosomes
#     gr_clipped = pr.gf.genome_bounds(gr_reduced, chromsize_gr, clip=True)
    

#     merged_ranges = gr.join(gr_clipped, suffix="_2")
    

#     merged_df = merged_ranges.sort(strand=False).df
#     merged_bg = merged_df.groupby(['Chromosome', 'Start_2', 'End_2']).agg({'values':['mean'],})
#     merged_bg.reset_index(inplace=True)

#     return merged_bg

#%%
# def interval_inclusion(bg_class, chrom_sizes):
#     '''
#     Pyranges overlap seems to do (1,5] while gr.merge is (1,5) 
#     Check if intervals are (0,10] or (1,9) format and make
#     all intervals in (1,9) non-overlapping format
#     Make sure coordinates don't run past chromosome end

#     Args:
#         bg_class (class): Bedgraph class
#         chrom_sizes (Ordered_dict): Chromosome sizes ordered dictionary from assembly.py

#     Note:
#         bedgraph should be 0-based coordinate system
#         can use genome_bounds for chromosome boundary

#         use slack=-1 to not merge bookended coordinates (x,y]
#     '''

#     overlap = set(bg_class.start).intersection(bg_class.end)
#     # only do comparison for fisrt 10000
#     if len(overlap) > 0:
#         print('Overlaping start end coordinates present')
#         # end_coord = [int(element)-1 for element in bg_class.end]
#         max_end = chrom_sizes.get(bg_class.chrom[0], 'missing')
#         current_chrom = bg_class.chrom[0]
#         for i in tqdm.tqdm(range(bg_class.length)):
#             chrom = bg_class.chrom[i]
#             if current_chrom != chrom:
#                 max_end = chrom_sizes.get(bg_class.chrom[i], 'missing')
#                 current_chrom = bg_class.chrom[i]

#             end = bg_class.end[i]
#             start = bg_class.start[i]
#             #modify end coordinates to end of chromosome
#             if max_end == 'missing':
#                 raise Exception('Chromosome notation incorrect')
#             elif end > max_end:
#                 bg_class.end[i] = int(max_end)
#             bg_class.start[i] = start+1
#             converted = True
#     else:
#         converted = False

#     return converted


 