import requests
import argparse
from collections import defaultdict
import re

def main():
    '''
    At the time of writing this the "full" format was broken so had 
    to use fasta and summary together

    # Repeat name (name, repName), 
    # Class (repeat_type_name, repClass), 
    # Family (repeat_subtype_name, repFamily), 
    # Other information (transcript name, accession)
    '''

    args = parse_arguments()

    if args.species == 'hs':
        clade = "Homo sapiens"
    elif args.species == 'mm':
        clade = "Mus musculus"
    else:
        print('Only mm or hs possible')
        sys.exit()

    dfam_v = dfam_version()

    clade_str = re.sub(' ', '_', clade)
    out_fa = clade_str + '_dfam_' + dfam_v + '.fasta'
    out_path = args.output + '/' + out_fa

    try:
        full = dfam_request("full", clade)
        process_full(full, out_path)
    except KeyError:
        summary = dfam_request("summary", clade)
        fa = dfam_request("fasta", clade)
        process_summary_fasta(summary, fa, out_path)

    print('Done')


def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Create a consensus fasta file from dfam')
    parser.add_argument('-o', '--output',
                        action = "store",
                        required=True,
                        help = "Output directory")
    parser.add_argument('-s', '--species',
                        type = str,
                        action = 'store',
                        required=True,
                        help = "mm or hs")

    return parser.parse_args()


def dfam_version():
    '''
    Get Dfam version being used to generate files
    '''
    
    url = "https://dfam.org/api/version"
    response = requests.get(url)
    version = '.'.join(list(response.json().values()))

    return version




def dfam_request(format, clade):
    '''Retreave query from the Dfam api
    
    Args:
        format(str): "fasta", "full", "summary"
        clade(str): "Mus musculus", "Homo sapiens"
    '''
    url = "https://dfam.org/api/families"
    params = {
        # The summary format is metadata-only and does not include
        # full details such as the consensus sequence and citations
        "format": format,

        # Only retrieve the first 10 results in this query
        # "limit": "10",

        # Search in Caenorhabditis elegans (worm)
        "clade": clade,

        # Include families from ancestor and descendant taxa in the results
        "clade_relatives": "both",
    }
    response = requests.get(url, params=params)
    if format == 'fasta':
        results = []
        for line in response.iter_lines():
            results.append(line)
    else:
        results = response.json()["results"]

    return results
    


def fasta_parse(fp):
    '''Parse fasta files
    Allows parsing of multiline fasta files

    Args:
        fp (str): Open fasta file
    '''

    record_out = 0
    record_count = 0
    name, seq = [''] * 2

    for line in fp:
        
        try:
            rec = line.decode('utf-8').rstrip()
        except AttributeError:
            rec = line.rstrip()

        if rec.startswith('>'):
            if record_count - record_out > 0:
                 yield name, seq
                 name, seq = [''] * 2

            name = rec
            record_count += 1
        elif rec.startswith('@'):
            raise Exception('Input starts with @ suggesting this is a fastq no fasta')
        else:
            seq += rec
    yield name, seq


def match_rmsk_names(name):
    '''
    Removes '5end', '3end', 'orf2' from Dfam names to match rmsk
    '''
    remove_ends = ['_5end', '_3end', '_orf2']

    for i in remove_ends: 
        if i in name:
            return(re.sub(i,'', name), i)
    
    return(name, '')


def process_full(results, fasta_out):
    '''
    Write out fasta with name rtn rsn and consensus sequence of repeats

    Args:
        results(list): dfam_request with full setting
        fasta_out(str): out path of consensus sequence fasta
    '''

    full_fields = ['accession', 'name', 'repeat_type_name', 'repeat_subtype_name', 
                  'consensus_sequence']

    #remove duplicate entries
    entries = defaultdict()
    dup = 0

    with open(fasta_out, 'w') as fa_out:
        for record in results:
            accession, name, rtname, rsname, con_seq = [record.get(i, 'NA') for i in full_fields]
            name_rmsk, end = match_rmsk_names(name)
            rec_name = '>' + name_rmsk + ',' + str(rtname) + ',' + str(rsname) + ',' + accession + end
            if rec_name in entries.keys():
                if con_seq == entries.get(rec_name):
                    dup += 1
            else:
                entries[rec_name] = con_seq
                fa_out.write(rec_name + '\n')
                seq_len = len(con_seq)
                written_out = 0
                while seq_len > 0:
                    fa_out.write(con_seq[written_out:60+written_out] + '\n')
                    written_out += 60
                    seq_len -= 60
    print('Duplicate sequences ignored:', dup)


def process_summary_fasta(results, fasta, fasta_out):
    '''
    Using summary and fasta results from dfam api, make lookup table and
    create a full fasta

    Args:
        results(list): dfam_request with summary setting
        fasta(str): dfam_request with fasta setting
        fasta_out(str): out path of consensus sequence fasta
    '''
    name_lookup = defaultdict()

    summary_fields = ['accession', 'name', 'repeat_type_name', 'repeat_subtype_name']

    for record in results:
        accession, name, rtn, rsn = [record.get(i, 'NA') for i in summary_fields]
        name_rmsk, end = match_rmsk_names(name)
        name_lookup[accession] = ','.join([name_rmsk, rtn, rsn, accession+end])

    with open(fasta_out, 'w') as fa_out:
        for rec, seq in fasta_parse(fasta):
            accs = rec.strip('>').split(' ')[0]
            fa_out.write('>' + name_lookup.get(accs) + '\n')
            seq_len = len(seq)
            written_out = 0
            while seq_len > 0:
                fa_out.write(seq[written_out:60+written_out] + '\n')
                written_out += 60
                seq_len -= 60



if __name__ == "__main__":
    main()


# #make fasta file with consensus_sequence
# with open('/mnt/data/mm10_repeats_consensus.fasta', 'w') as fa_out:
#     for record in results:
#         name, rtname, rsname, con_seq = [record[i] for i in ['name',  'repeat_type_name', 'repeat_subtype_name', 'consensus_sequence']]
#         fa_out.write('>' + name + '|' + str(rtname) + '|' + str(rsname) + '\n')
#         seq_len = len(con_seq)
#         written_out = 0
#         while seq_len > 0:
#             fa_out.write(con_seq[written_out:80+written_out] + '\n')
#             written_out += 80
#             seq_len -= 80
# # Prints "Vingi-2_CE" at the time of this writing
# # print(results[2]["name"])



# with open('/mnt/data/mm10_repeats_names.txt', 'w') as txt_out:
#     for record in results:
#         name, repeat_type_name, repeat_subtype_name = [record[i] for i in ['name', 'repeat_type_name', 'repeat_subtype_name']]
#         txt_out.write(str(name) + '\t' + str(repeat_type_name) + '\t' + str(repeat_subtype_name) + '\n')

# for record in results:
#     if 'RLTR17B_Mm' in record.values():
#         print(record)

#find distance between repeats



# params_hs = {
#     # The summary format is metadata-only and does not include
#     # full details such as the consensus sequence and citations
#     "format": "full",

#     # Only retrieve the first 10 results in this query
#     # "limit": "10",

#     # Search in Caenorhabditis elegans (worm)
#     "clade": "Homo sapiens",

#     # Include families from ancestor and descendant taxa in the results
#     "clade_relatives": "both",
# }


# response_hs = requests.get(url, params=params_hs)
# results_hs = response_hs.json()["results"]

# #make fasta file with consensus_sequence
# with open('/mnt/data/hs_repeats_consensus.fasta', 'w') as fa_out:
#     for record in results_hs:
#         name, rtname, rsname, con_seq = [record[i] for i in ['name',  'repeat_type_name', 'repeat_subtype_name', 'consensus_sequence']]
#         fa_out.write('>' + name + '|' + str(rtname) + '|' + str(rsname) + '\n')
#         seq_len = len(con_seq)
#         written_out = 0
#         while seq_len > 0:
#             fa_out.write(con_seq[written_out:80+written_out] + '\n')
#             written_out += 80
#             seq_len -= 80

