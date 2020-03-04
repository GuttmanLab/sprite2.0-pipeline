
import gzip
import re
from collections import Counter, defaultdict
import pandas as pd
import argparse
import seaborn as sns
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import os
from functools import reduce
import fnmatch

'''
1. Quantify the proportions of the 12 different barcodes added at each round.
   Barcode 12 was left empty to determine cross-ligation
2. Count number of reads with all 8 identical barcodes to estimate efficiency of
   capturing full complexes.
3. Count barcodes at each ligation round, with 96 unique barcodes can determine if
   incorrect sticky ends ligated together
'''


# fastqgzfile = '/mnt/data/20190621_PB_PC_SQ/assembled/PB-PC_S2_L001_R1_001_val_1.assembled_rv_bID.fq.gz'
# out_table = ''


def parse_args():

    parser = argparse.ArgumentParser(description='Quantify proportions of the 12 barcodes in each of the 8 barcoding rounds')
    parser.add_argument('--r1', dest='read_1', type=str, required=True,
                        help='Fastq read 1')
    parser.add_argument('--to', dest='table', type=str, required=True,
                        help='Count table out')
    parser.add_argument('--single', dest='single', action='store_true',
                        help='Single plate')
    # parser.add_argument('--sample', dest='sample', type=str, required=False, default=False,
    #                     help='R1 is used as sample barcode')

    opts = parser.parse_args()

    return opts


def main():

    #argparse
    opts = parse_args()

    pattern = re.compile('\[([a-zA-Z0-9_\-]+)\]')
    cluster_sizes = {'1000x':Counter(), '100x':Counter(), '10x':Counter(), '1x':Counter(), 'empty':Counter()}
    # cluster_sizes = defaultdict(int)
    round_1 = {1:Counter(),2:Counter(),3:Counter(),4:Counter(),5:Counter(),6:Counter(),7:Counter(),8:Counter()}
    round_2 = {1:Counter(),2:Counter(),3:Counter(),4:Counter(),5:Counter(),6:Counter(),7:Counter(),8:Counter()}
    round_3 = {1:Counter(),2:Counter(),3:Counter(),4:Counter(),5:Counter(),6:Counter(),7:Counter(),8:Counter()}
    sticky_ends= defaultdict(lambda: defaultdict(int))
    sticky_ends_full = defaultdict(lambda: defaultdict(int))

    barcodes_lst = []
    total = 0
    skipped = 0

    with gzip.open(opts.read_1, "rt") as f:
    # with gzip.open(fastqgzfile, "rt") as f:
        for line in f:
            barcodes = pattern.findall(line)

        # if total <10:
            #Sticky-end quantification - made possible by using 96 unique barcodes
                # ['TermStag_bot_6', 'R7_bot_A2', 'R6_bot_B10', 'R5_bot_A7', 'R4_bot_A9', 'R3_bot_A8', 'NOT_FOUND', 'NOT_FOUND']
            for r in range(len(barcodes)):
                sticky_ends['position_' + str(r)][barcodes[r].split('_')[0]] += 1

            if 'NOT_FOUND' not in barcodes:

                barcodes_all = '|'.join(barcodes)
                barcodes_lst.append(barcodes_all)

                if opts.single:
                    cluster_size = barcodes[7].split('_')[-1].strip('AB')
                    # cluster_sizes[cluster_size][barcodes_all] += 1
                    if cluster_size in ['1','2','3','4']:
                        cluster_sizes['1000x'][barcodes_all] += 1
                    elif cluster_size in ['5','6','7']:
                        cluster_sizes['100x'][barcodes_all] += 1
                    elif cluster_size in ['8','9']:
                        cluster_sizes['10x'][barcodes_all] += 1
                    elif cluster_size in ['10','11']:
                        cluster_sizes['1x'][barcodes_all] += 1
                    elif cluster_size == '12':
                        cluster_sizes['empty'][barcodes_all] += 1
                    else:
                        print('no found', cluster_size)
                else:
                    #create a string to group barcodes by
                    r1 = [barcodes[-1]]

                    if len(fnmatch.filter(r1, 'R1-*_bot_A[1-4]')) > 0:
                        cluster_sizes['1000x'][barcodes_all] += 1
                    elif len(fnmatch.filter(r1, 'R1-*_bot_A[5-7]')) > 0:
                        cluster_sizes['100x'][barcodes_all] += 1
                    elif len(fnmatch.filter(r1, 'R1-*_bot_A[8-9]')) > 0:
                        cluster_sizes['10x'][barcodes_all] += 1
                    elif len(fnmatch.filter(r1, 'R1-*_bot_A1[0-1]')) > 0:
                        cluster_sizes['1x'][barcodes_all] += 1
                    elif len(fnmatch.filter(r1, 'R1-*_bot_A12')) > 0:
                        cluster_sizes['empty'][barcodes_all] += 1
                    else:
                        print('unassigned barcodes', r1)



                #go through all barcodes to count
                for r in range(len(barcodes)):
                    if opts.single:
                        round_1[r+1][barcodes[r].split('_')[-1].strip('AB')] += 1
                    else:
                        plate = barcodes[r].split('-')[-1][0]
                        if plate is '1':
                            round_1[r+1][barcodes[r].split('_')[-1].strip('AB')] += 1
                        elif plate is '2':
                            round_2[r+1][barcodes[r].split('_')[-1].strip('AB')] += 1
                        elif plate is '3':
                            round_3[r+1][barcodes[r].split('_')[-1].strip('AB')] += 1
                        else:
                            print('plate not found', plate)

                    # round[r+1][barcodes[r].split('_')[-1].strip('AB')] += 1
                    sticky_ends_full['position_' + str(r)][barcodes[r].split('_')[0]] += 1
            else:
                skipped += 1

            next(f)
            next(f)
            next(f)
            total += 1
    # else:
    #     break


    if opts.single:
        final_df_1 = pd.DataFrame.from_dict(round_1, orient='index')
        col_order = sorted(final_df_1.columns.astype('int').to_list())
        col_order = list(map(str, col_order))
        final_df_1 = final_df_1[col_order]
        final_df_1.to_csv(os.path.splitext(opts.table)[0] + '_all_plates.csv', sep='\t')
    else:
        #write out counts table
        final_df_1 = pd.DataFrame.from_dict(round_1, orient='index')
        final_df_2 = pd.DataFrame.from_dict(round_2, orient='index')
        final_df_3 = pd.DataFrame.from_dict(round_3, orient='index')



        dfs = [final_df_1, final_df_2, final_df_3]
        final_df = reduce(lambda x, y: x.add(y, fill_value=0), dfs)
        # final_df.index = final_df.index.astype(int)
        # final_df.sort_index(axis=1, inplace=True)
        # col_order = sorted(final_df.columns.astype('int').to_list())
        # col_order = list(map(str, col_order))
        # final_df = final_df[col_order]
        # final_df.to_csv(opts.table, sep='\t')


        write_table(final_df_1, os.path.splitext(opts.table)[0] + '_plate1.csv')
        write_table(final_df_2, os.path.splitext(opts.table)[0] + '_plate2.csv')
        write_table(final_df_3, os.path.splitext(opts.table)[0] + '_plate3.csv')

    #write out sticky end table
    sticky_table = pd.DataFrame.from_dict(sticky_ends, orient='index')
    col_order = sorted(sticky_table.columns.to_list())
    col_order = col_order[1:] + [col_order[0]]
    sticky_table = sticky_table[col_order]
    sticky_table.to_csv(os.path.splitext(opts.table)[0] + '_sticky_lig.csv', sep='\t')

    #write out sticky end table full
    sticky_table_full = pd.DataFrame.from_dict(sticky_ends_full, orient='index')
    col_order_full = sorted(sticky_table_full.columns.to_list())
    sticky_table_full = sticky_table_full[col_order_full]
    sticky_table_full.to_csv(os.path.splitext(opts.table)[0] + '_sticky_lig_full.csv', sep='\t')

    print('Total: ' + str(total))
    print('Skipped: ' + str(skipped))

    #make violin plot of cluster sizes
    # with PdfPages(opts.table.split('.')[0] + '_violinplot.pdf') as pdf_out:
    pd_1000x = pd.DataFrame({'cluster':np.array(['1000x']*len(np.array(list(cluster_sizes['1000x'].values())))),
                            'size':np.array(list(cluster_sizes['1000x'].values()))})

    pd_100x = pd.DataFrame({'cluster':np.array(['100x']*len(np.array(list(cluster_sizes['100x'].values())))),
                            'size':np.array(list(cluster_sizes['100x'].values()))})

    pd_10x = pd.DataFrame({'cluster':np.array(['10x']*len(np.array(list(cluster_sizes['10x'].values())))),
                           'size':np.array(list(cluster_sizes['10x'].values()))})

    pd_1x = pd.DataFrame({'cluster':np.array(['1x']*len(np.array(list(cluster_sizes['1x'].values())))),
                            'size':np.array(list(cluster_sizes['1x'].values()))})

    pd_empty = pd.DataFrame({'cluster':np.array(['empty']*len(np.array(list(cluster_sizes['empty'].values())))),
                            'size':np.array(list(cluster_sizes['empty'].values()))})

    plot_df = pd.concat([pd_1000x, pd_100x, pd_10x, pd_1x, pd_empty])


        # sns.set(style="whitegrid")
        # violin_plot = sns.violinplot(data=plot_df, x='cluster', y='size', inner='boxplot')
        # fig = violin_plot.get_figure()
        # pdf_out.savefig(fig)

    plot_df.to_csv(os.path.splitext(opts.table)[0] + '_hist.txt', sep='\t')

    with open(os.path.splitext(opts.table)[0] + '_barcodes.txt', 'w') as out_br:
        for item in barcodes_lst:
            out_br.write(item + '\n')

def write_table(df, out_path):
    '''Write out a pandas dataframe to csv

    Args:
        df (DataFrame): pandas dataframe to write out
        out_path (str): out path
    '''
    col_order = sorted(df.columns.astype('int').to_list())
    col_order = list(map(str, col_order))
    df = df[col_order]
    df.to_csv(out_path, sep='\t')




if __name__ == "__main__":
    main()
