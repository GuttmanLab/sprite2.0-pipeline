
from babrahamlinkon.general import file_open, fastq_parse
import numpy as np
import pandas as pd
import seaborn as sns
import regex
import re
from colorama import Fore, Style
from collections import defaultdict
import argparse

# in_fq_path = '/mnt/data/SPRITEzero_complexes/trimmed/assembled/25000_S1_L001_R1_001_val_1_assembled_rv_bID.fastq.gz'
# in_fq_path_new = '/mnt/data/20190509_newadaptsprite0/workup/assembled/100_S1_L001_R1_001_val_1_assembled_rv.fastq_bID.fq.gz'

# seq_lengths = []
# #examine read length
# with file_open(in_fq_path_new) as in_fq:
#     for qname, seq, thrd, qual in fastq_parse(in_fq):
#          seq_lengths.append(len(seq))
#
#
# x = pd.Series(seq_lengths, name='Seq lengths')
# ax = sns.distplot(x)

def parse_args():

    parser = argparse.ArgumentParser(description='Examine sticky ends within sequences')
    parser.add_argument('--r2', dest='read_2', type=str, required=True,
                        help='Fastq read 2 - looks for bottom sticky ends')
    parser.add_argument('--max', dest='max', type=int, required=False, default=0,
                        help='Either the number of barcodes (examine sequences with only 4 barcodes) \
                        or max length of sequences to examine. (0)')
    parser.add_argument('--min', dest='min', type=int, required=False,
                        default=0, help='Min length of sequences to examine. (0)')
    parser.add_argument('--tags', dest='tags', type=int, required=False,
                        default=8, help='Number of barcoding rounds/tags (8)')
    parser.add_argument('--printn', dest='printn', type=int, required=False,
                        default=100, help='Number of annotated sequences to print (100)')
    opts = parser.parse_args()

    return opts


def main():

    opts = parse_args()

    seq_dict = get_seq_to_anno(opts.read_2, min=opts.min, max=opts.max, number_of_barcodes=opts.tags)
    print_anno_seq(seq_dict, opts.printn)



def get_seq_to_anno(fastq_path, min=0, max=0, number_of_barcodes=8):
    '''Get reads of interest (i.e. missing barcodes, short length)
    '''
    #get reads with 4 and 6 Barcodes
    # four_bc = defaultdict(list)
    # six_bc = defaultdict(list)
    # all_seqs = defaultdict(list)
    # right_len = defaultdict(list)

    seq_dict = defaultdict(list)

    pattern = re.compile('\[([a-zA-Z0-9_\-]+)\]')

    with file_open(fastq_path) as in_fq:
        for qname, seq, thrd, qual in fastq_parse(in_fq):
            barcodes = pattern.findall(qname)
            num_mis_bc = barcodes.count('NOT_FOUND')

            if max > number_of_barcodes: #if more than this is sequence length, else its num barcodes
                if len(seq) < max and len(seq) > min:
                    seq_dict['|'.join(barcodes)].append(seq)

            elif num_mis_bc == max:
                seq_dict['|'.join(barcodes)].append(seq)


    print('Captured sequences:', sum(len(v) for v in seq_dict.values()))

    return seq_dict

# 1648911/2345224
#0.7030931800118028 correct length

def print_anno_seq(seq_dict, num_print=100):
    '''Print annotation of sequences in dictionary

    Args:
        seq_dict (dict): A dictionary of barcodes (key) and sequences (values)
        num_print (int): number of sequences to print
    '''

    count = 0
    for k, v in seq_dict.items():
        anno_seqs = []
        if count < num_print:
            for seq in v:
                if count < num_print:
                    anno_seqs.append(sticky2ansi(seq))
                    count += 1
                else:
                    break
        else:
            break
        print(k + '\n' + '\n'.join(anno_seqs))



# bc_exl = pd.ExcelFile('/mnt/data/Google Drive/Labbook/SPRITE/SPRITEzero barcodes design/20190118_R1-R8_8mers.xlsx')
# bc_exl = pd.ExcelFile('/mnt/data/Google Drive/Labbook/SPRITE/SPRITEzero barcodes design/20190502_R1-R8_8mers_spritezero.xlsx')
#
# print(bc_exl.sheet_names)
# bc_df = bc_exl.parse('Plate')
#
# bc_dict = pd.DataFrame.to_dict(bc_df.iloc[:,0:3], orient='records')
#
# barcodes_dict = defaultdict(lambda: {'sticky_end':'', 'barcode':''})
# for rec in bc_dict:
#     name = rec['Unnamed: 0']
#     sticky_end = rec['sticky end']
#     barcode = rec['barcode']
#
#     barcodes_dict[name]['sticky_end'] = sticky_end
#     barcodes_dict[name]['barcode'] = barcode
#
# config = '/mnt/data/20190509_newadaptsprite0/workup/assembled/config.txt'
#
# barcodes_dct = defaultdict()
# with open(config, 'r') as cf:
#     for _ in range(4):
#         next(cf)
#     for line in cf:
#         tag, name, barcode, errors = line.rstrip().split('\t')
#         barcodes_dct[name] = barcode
#
# #
# # from colorama import Fore, Style
# # col = 'GREEN'
# # print(f'{Fore.WHITE}test {Fore.GREEN}color{Style.RESET_ALL}!')
# #
# # Fore.GREEN
#
# anno_keys = list(right_len.keys())[0].split('|')
# str_to_annotate = list(right_len.values())[0][0]
#
# bc_to_anno = [barcodes_dct.get(bc, 'NOT_FOUND') for bc in anno_keys]


def str2ansi(barcodes, str_to_annotate):
    '''Add ANSI codes into string based on location of barcodes

    '''
    bc_colors = [Fore.MAGENTA, Fore.RED, Fore.YELLOW, Fore.BLACK, Fore.BLUE,
                 Fore.CYAN, Fore.GREEN, Fore.LIGHTBLACK_EX]

    barcode = 0
    str_out = str_to_annotate
    for bc in barcodes:
        assert barcode < 8, 'More than 8 barcodes present'
        #skip annotation if bc is not found
        if bc is 'NOT_FOUND':
            barcode += 1
            continue
        else:
            try:
                bc_found = regex.search(bc + "{e<=1}", str_out) # means allow up to 1 error
                start = bc_found.start()
                end = bc_found.end()
                assert end - start == len(bc), 'barcode length incorrect'
                str_out = str_out[:start] + bc_colors[barcode] + str_out[start:end] + Style.RESET_ALL + str_out[end:]
                barcode += 1
            except:
                continue
    # assert len(str_out) == len(str_to_annotate), 'Out string is a different length then input string!'

    return str_out


# print(str2ansi(bc_to_anno, str_to_annotate))
# regex.search(bc_to_anno[1] + "{e<=1}", str_to_annotate)




def sticky2ansi(seq):
    '''Add ANSI codes into string based on location of sticky ends
    '''
    bc_colors = [Fore.MAGENTA, Fore.RED, Fore.YELLOW, Fore.BLACK, Fore.BLUE,
                 Fore.CYAN, Fore.GREEN, Fore.LIGHTBLACK_EX]

    # [TermStag_bot_4][R7_bot_A7][R6_bot_B5][R5_bot_A1][R4_bot_A5][R3_bot_A5][R2_bot_B6][R1_A9_10x_9]
    # text = """TCGTTCGACGGCATACCCAGCGGAGAGCGTTGCCGACCTTTGACGTCGGCAAGCGCTGATAGACCGGCTATCTGCTGACGCGAGCAACAGCCCGGTTCCACGAGAGGGCGGTGATGACTTGGTTCCCAGGGGTCGGC"""

    # regex.search("(?e)(dog){e<=1}", "cat and dog")[1] returns "dog"
    # (without a leading space) because the fuzzy search matches " dog" with 1 error,
    # which is within the limit, and the (?e) then it attempts a better fit.
    regex_se = regex.compile(r"(?e)(TGACTTG){e<=1}|(?e)(ACGAGAG){e<=1}|(?e)(CAACAGC){e<=1}|(?e)(ATCTGCT){e<=1}|(?e)(GCTGATA){e<=1}|(?e)(TTGACGT){e<=1}|(?e)(GAGCGTT){e<=1}|(?e)(GGCATAC){e<=1}", regex.I)

    i = 0; output = ''
    for m in regex_se.finditer(seq):
        output += ''.join([seq[i:m.start()],
                           bc_colors[m.lastindex-1],
                           seq[m.start():m.end()],
                           Style.RESET_ALL,])
        i = m.end()

    return ''.join([output, seq[m.end():]])



if __name__ == "__main__":
    main()
