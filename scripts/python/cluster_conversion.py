import pandas as pd
import numpy as np
import sys
import re
import tqdm

# Cluster file
clusterfile = sys.argv[1]

# Output file
outfile = sys.argv[2]

# Minsize
minsize = int(sys.argv[3])

# Maxsize
maxsize = int(sys.argv[4])

# Conversion type
contype = sys.argv[5]


###########################################################################
# Conversion Functions
###########################################################################

def convert_clusters_dpm(clusters, output, min_size, max_size):
    '''
    Convert clusters containing DPM reads
    '''
    store = pd.HDFStore(output, 'w', complevel=0, complib='blosc') #lzo
    df = parse_dpm(clusters)
    barcodes = parse_barcodes(clusters)
    barcodes = barcodes.query('DPM!=0')
    # Filter for size
    barcodes  = barcodes.query('Size>=@min_size and Size<=@max_size')
    df = df.query('Size>=@min_size and Size<=@max_size')
    # Reassign the barcodes
    barcodes = barcodes.reset_index(drop=True)
    barcode_dict = barcodes.reset_index().set_index('Name')['index'].to_dict()
    barcodes = barcodes.reset_index()
    barcodes.rename(columns = {'Name':'Barcode', 'index':'Name'}, inplace = True)
    print(barcodes.head(100))
    df = simplify_barcodes(barcode_dict,df)
    print(df.head(100))
    # Split the dpm reads
    dpm_only = df.query('Size==DPMSize')
    dpm_rpm = df.query('Size!=DPMSize')
    # Add to the store
    store.append('DPM_only', dpm_only, format='table', append=True, data_columns=True, index=False)
    store.append('DPM_RPM', dpm_rpm, format='table', append=True, data_columns=True, index=False)    
    store.append('Names', barcodes, format='table', append=True, data_columns=True, index=False)
    store.create_table_index('DPM_only', columns=list(dpm_only.columns), optlevel=9, kind='full')
    store.create_table_index('DPM_RPM', columns=list(dpm_rpm.columns), optlevel=9, kind='full')
    store.create_table_index('Names', columns=list(barcodes.columns), optlevel=9, kind='full')
    store.close()


def convert_clusters_rpm(clusters, output, min_size, max_size):
    '''
    Convert clusters containing only RNA reads
    '''
    store = pd.HDFStore(output, 'w', complevel=0, complib='blosc') #lzo
    df = parse_rpm(clusters)
    repeats = parse_repeats(clusters)
    barcodes = parse_barcodes(clusters)
    # Filter for rna barcodes without dna
    barcodes = barcodes.query('RPM!=0 and DPM==0')
    df = df.query('DPMSize==0')
    repeats = repeats.query('DPMSize==0')
    # Filter for size
    barcodes  = barcodes.query('Size>=@min_size and Size<=@max_size')
    df = df.query('Size>=@min_size and Size<=@max_size')
    repeats = repeats.query('Size>=@min_size and Size<=@max_size')
    # Reassign the barcodes
    barcodes = barcodes.reset_index(drop=True)
    barcode_dict = barcodes.reset_index().set_index('Name')['index'].to_dict()
    barcodes = barcodes.reset_index()
    barcodes.rename(columns = {'Name':'Barcode', 'index':'Name'}, inplace = True)
    print(barcodes.head(100))
    df = simplify_barcodes(barcode_dict,df)
    repeats = simplify_barcodes(barcode_dict, repeats, True)
    print(repeats.head(100))
    # Parse features of RPM
    df = split_feature(df)
    print(df.head(100))
    # Add to the store
    store.append('RPM', df, format='table', append=True, data_columns=True, index=False)
    store.append('RepeatRPM', repeats, format='table', append=True, data_columns=True, index=False)    
    store.append('Names', barcodes, format='table', append=True, data_columns=True, index=False)
    store.create_table_index('RPM', columns=list(df.columns), optlevel=9, kind='full')
    store.create_table_index('RepeatRPM', columns=list(repeats.columns), optlevel=9, kind='full') 
    store.create_table_index('Names', columns=list(barcodes.columns), optlevel=9, kind='full')
    store.close()


def convert_clusters_rpm_dpm(clusters, output, min_size, max_size):
    '''
    Convert clusters containing both RNA and DNA reads
    '''
    store = pd.HDFStore(output, 'w', complevel=0, complib='blosc') #lzo
    rpm = parse_rpm(clusters)
    dpm = parse_dpm(clusters)
    repeats = parse_repeats(clusters)
    barcodes = parse_barcodes(clusters)
    # Filter for rna-dna
    barcodes = barcodes.query('RPM!=0 and DPM!=0')
    rpm = rpm.query('DPMSize!=0')
    repeats = repeats.query('DPMSize!=0')
    dpm = dpm.query('DPMSize!=Size')
    # Filter for size
    barcodes  = barcodes.query('Size>=@min_size and Size<=@max_size')
    dpm = dpm.query('Size>=@min_size and Size<=@max_size')
    rpm = rpm.query('Size>=@min_size and Size<=@max_size') 
    repeats = repeats.query('Size>=@min_size and Size<=@max_size')
    # Reassign the barcodes
    barcodes = barcodes.reset_index(drop=True)
    barcode_dict = barcodes.reset_index().set_index('Name')['index'].to_dict()
    barcodes = barcodes.reset_index()
    barcodes.rename(columns = {'Name':'Barcode', 'index':'Name'}, inplace = True)
    print(barcodes.head(100))
    dpm = simplify_barcodes(barcode_dict,dpm)
    rpm = simplify_barcodes(barcode_dict,rpm)
    repeats = simplify_barcodes(barcode_dict, repeats, True)
    print(dpm.head(100))
    print(repeats.head(100))
    # Annotate RNAs
    rpm = split_feature(rpm) 
    print(rpm.head(100))
    # Add to the store
    store.append('RPM', rpm, format='table', append=True, data_columns=True, index=False)
    store.append('DPM', dpm, format='table', append=True, data_columns=True, index=False)
    store.append('RepeatRPM', repeats, format='table', append=True, data_columns=True, index=False)    
    store.append('Names', barcodes, format='table', append=True, data_columns=True, index=False)
    store.create_table_index('RPM', columns=list(rpm.columns), optlevel=9, kind='full')
    store.create_table_index('DPM', columns=list(dpm.columns), optlevel=9, kind='full')
    store.create_table_index('RepeatRPM', columns=list(repeats.columns), optlevel=9, kind='full')
    store.create_table_index('Names', columns=list(barcodes.columns), optlevel=9, kind='full')
    store.close()

##########################################################################
#Parsing Functions
##########################################################################

def parse_barcodes(clusters):
    '''
    Count number of DPM and RPM reads associated with each cluster
    '''
    barcodes, rpm, dpm, size = list(), list(), list(), list()
    with open(clusters, "r") as c:
        for line in tqdm.tqdm(c):
            barcode, *reads = line.rstrip('\n').split('\t')
            barcodes.append(barcode)
            total_size = len(reads)
            dpm_reads = [r for r in reads if r.startswith("DPM")]
            dpm_size = len(dpm_reads)
            size.append(total_size)
            dpm.append(dpm_size)
            rpm.append(total_size - dpm_size)
    df = pd.DataFrame({'Name': pd.Series(barcodes, dtype=str),
                       'RPM': pd.Series(rpm, dtype=int),
                       'DPM': pd.Series(dpm, dtype=int),
                       'Size': pd.Series(size, dtype=int)})
    return df


def parse_dpm(clusters):
    '''
    Parse DPM reads from cluster into a dataframe
    '''
    barcodes, chroms, starts, ends, strands, size, dpm_sizes = list(), list(), list(), list(), list(), list(), list()
    pattern = re.compile('([a-zA-Z0-9]+)\[(.+)\]_(.+):([0-9]+)\-([0-9]+)')
    with open(clusters, "r") as c:
        for line in tqdm.tqdm(c):
            barcode, *reads = line.rstrip('\n').split('\t')
            total_size = len(reads)
            dpm_reads = [r for r in reads if r.startswith('DPM')]
            dpm_size = len(dpm_reads)
            for read in dpm_reads:
                match = pattern.search(read)
                read_type, feature, chrom, start, end = match.groups()
                size.append(total_size)
                chroms.append(chrom)
                strands.append(feature)
                starts.append(int(start))
                ends.append(int(end))
                dpm_sizes.append(dpm_size)
                barcodes.append(barcode)
    df = pd.DataFrame({'Chromosome': pd.Series(chroms, dtype=str),
                       'Start': pd.Series(starts, dtype=int),
                       'End': pd.Series(ends, dtype=int),
                       'Name': pd.Series(barcodes, dtype=str),
                       'Strand': pd.Series(strands, dtype=str),
                       'Size': pd.Series(size, dtype=int),
                       'DPMSize': pd.Series(dpm_sizes, dtype=int)})
    return df


def parse_rpm(clusters):
    '''
    Parse RNA reads that aligned to the genome from cluster file
    '''
    barcodes, chroms, starts, ends, strands, features, size, dpm_sizes = list(), list(), list(), list(), list(), list(), list(), list()
    pattern = re.compile('([a-zA-Z0-9]+)\[(.+)\]_(.+):([0-9]+)\-([0-9]+)')
    with open(clusters, "r") as c:
        for line in tqdm.tqdm(c):
            barcode, *reads = line.rstrip('\n').split('\t')
            total_size = len(reads)
            rpm_reads = [r for r in reads if r.startswith("RPM")]
            dpm_size = total_size- len(rpm_reads)
            for read in rpm_reads:
                match = pattern.search(read)
                read_type, feature, chrom, start, end = match.groups()
                if chrom.startswith('chr'):
                    size.append(total_size)
                    chroms.append(chrom)
                    anno, strand = feature.rsplit(';', 1)
                    strands.append(strand)
                    starts.append(int(start))
                    ends.append(int(end))
                    dpm_sizes.append(dpm_size)
                    barcodes.append(barcode)
                    features.append(anno)
    df = pd.DataFrame({'Chromosome': pd.Series(chroms, dtype=str),
                       'Start': pd.Series(starts, dtype=int),
                       'End': pd.Series(ends, dtype=int),
                       'Name': pd.Series(barcodes, dtype=str),
                       'Strand': pd.Series(strands, dtype=str),
                       'Feature': pd.Series(features, dtype=str),
                       'Size': pd.Series(size, dtype=int),
                       'DPMSize': pd.Series(dpm_sizes, dtype=int)})
    return df


def parse_repeats(clusters):
    '''
    Parse RNA reads that aligned to custom genome from cluster file
    '''
    barcodes, repeats, size, dpm_sizes  = list(), list(), list(), list()
    pattern = re.compile('([a-zA-Z0-9]+)\[(.+)\]_(.+):([0-9]+)\-([0-9]+)')
    with open(clusters, "r") as c:
        for line in tqdm.tqdm(c):
            barcode, *reads = line.rstrip('\n').split('\t')
            total_size = len(reads)
            rpm = [r for r in reads if r.startswith("RPM")]
            dpm_size =total_size- len(rpm)
            for read in rpm:
                match = pattern.search(read)
                read_type, feature, chrom, start, end = match.groups()
                if not chrom.startswith('chr'):
                    size.append(total_size)
                    dpm_sizes.append(dpm_size)
                    barcodes.append(barcode)
                    repeats.append(chrom)
    df = pd.DataFrame({'Name': pd.Series(barcodes, dtype=str),
                       'Repeat': pd.Series(repeats, dtype=str),
                       'Size': pd.Series(size, dtype=int),
                       'DPMSize': pd.Series(dpm_sizes, dtype=int)})
    return df


def simplify_barcodes(barcode_dict, df, repeat=False):
    '''
    Convert the barcode names into simple integers
    '''
    df.loc[:,'Name'] = df['Name'].map(barcode_dict)
    df = df.astype({'Name':'int64'}, copy=False)
    if repeat == True:
       df.sort_values(by=['Name', 'Repeat'], inplace=True)
    else:
       df.sort_values(by=['Name', 'Chromosome'], inplace=True)
    return df


def parse_feature(feature_string):
    '''Split annotation feature into three columns'''
    feature_dict = {}
    for item in feature_string.split(';'):
        anno, tag = item.rsplit('.',1)
        if tag in feature_dict:
            feature_dict[tag] = 'AMB'
        else:
            feature_dict[tag] = anno
    return (feature_dict.get('exon', str(np.nan)),
            feature_dict.get('intron', str(np.nan)),
            feature_dict.get('repeat', str(np.nan)))

def split_feature(df):
   '''Split annotation feature into three columns'''
   df['Split'] = df['Feature'].apply(parse_feature)
   df[['Exon', 'Intron', 'Repeat']] =  pd.DataFrame(df['Split'].tolist(), index=df.index)
   del df['Split']
   del df['Feature']
   return df


##################################################################
#Function Calls
################################################################


if contype == 'DPM':
    convert_clusters_dpm(clusterfile, outfile, minsize, maxsize)
elif contype == 'RPM':
    convert_clusters_rpm(clusterfile, outfile, minsize, maxsize)
elif contype == 'DPM_RPM':
    convert_clusters_rpm_dpm(clusterfile, outfile, minsize, maxsize)

