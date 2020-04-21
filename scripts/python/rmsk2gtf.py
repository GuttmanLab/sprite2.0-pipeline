import gzip
from collections import defaultdict
from tqdm import tqdm

#%%
def main():
    rmsk2gtf('/mnt/data/genomes/GRCm38.p6/repeats/rmskOutCurrent.txt.gz',
             '/mnt/data/genomes/GRCm38.p6/repeats/mm10_rmskOutCurrent.gtf.gz')

    rmsk2gtf('/mnt/data/genomes/GRCh38/repeats/rmsk.txt.gz',
             '/mnt/data/genomes/GRCh38/repeats/hg38_rmsk.gtf.gz')


#%%
'''
rmsk format:

field   example SQL type    description
bin 585 smallint(5) unsigned    Indexing field to speed chromosome range queries.
swScore 463 int(10) unsigned    Smith Waterman alignment score
milliDiv    13  int(10) unsigned    Base mismatches in parts per thousand
milliDel    6   int(10) unsigned    Bases deleted in parts per thousand
milliIns    17  int(10) unsigned    Bases inserted in parts per thousand
genoName    chr1    varchar(255)    Genomic sequence name
genoStart   10000   int(10) unsigned    Start in genomic sequence
genoEnd 10468   int(10) unsigned    End in genomic sequence
genoLeft    -248945954  int(11) -#bases after match in genomic sequence
strand  +   char(1) Relative orientation + or -
repName (TAACCC)n   varchar(255)    Name of repeat
repClass    Simple_repeat   varchar(255)    Class of repeat
repFamily   Simple_repeat   varchar(255)    Family of repeat
repStart    1   int(11) Start (if strand is +) or -#bases after match (if strand is -) in repeat sequence
repEnd  471 int(11) End in repeat sequence
repLeft 0   int(11) -#bases after match (if strand is +) or start (if strand is -) in repeat sequence
id  1   char(1) First digit of id field in RepeatMasker .out file. Best ignored.

UCSC rmsk is 0-based coordinate system
GFF is 1-based coordinate system
gtf format:
1       mm10_rmsk       exon    3044870 3045464 4839    +       .       gene_id "RLTR1D"; transcript_id "RLTR1D"; family_id "ERV1"; class_id "LTR";
1       mm10_rmsk       exon    3014472 3014749 1604    +       .       gene_id "L1VL1"; transcript_id "L1VL1"; family_id "L1"; class_id "LINE"; all_id "L1VL1,L1VL1,L1,LINE"
'''

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
def gft_format_meta(meta_data):
    '''Format a dictionary into gtf metadata format
    
    Args:
        meta_data(dict): dictionary with metadata
    '''
    out_f = []
    for k, v in meta_data.items():
        out_f.append(str(k) + ' ' + '"' + str(v) + '"')
    out = ';'.join(out_f)
    return out

#%%
class genome_generic:
    '''A generic class to hold genomic information

    Methods:
    - add_meta(): Add a dictionary with additional metadata
    "{nameOfFeature:'myFeature'}"
    '''

    def __init__(self, chrom, start, end, strand):
        self._chrom = chrom
        self._start = start
        self._end = end
        self._stand = strand
        self._meta = defaultdict(set)

    def add_meta(self, metadata):
        if isinstance(metadata, dict):
            for k, v in metadata.items():
                #each feature name can only have a single value
                self._meta[k] = v

    def to_gtf(self, name, gtype, original_data_0_based=True):
        '''
        conver to 1-based if original data is 0-based
        '''
        if original_data_0_based:
            start = int(self._start) + 1
        else:
            start = int(self._start)
        
        out_format = '\t'.join([self._chrom, name, gtype, str(start), self._end, 
        self._meta.get('score', '0'), self._stand, '.', gft_format_meta(self._meta), '\n'])
        return out_format


#%%
def rmsk2gtf(rmsk_path, gtf_out):
    '''Convert a ucsc repeatmask to a gtf format

    Args:
        rmsk_path(str): path to repeatmask file
        gtf_out(str): output path of gtf file
    '''
    count = 0
    keep_fields = [1, 2, 5, 6, 7, 9, 10, 11, 12]
    transcript_id_unique = defaultdict(int)
    with file_open(rmsk_path) as rmsk, \
        gzip.open(gtf_out, 'wt') as out_gtf:
        for line in tqdm(rmsk):
            # if count < 10:
            score, milliDiv, genoName, genoStart, genoEnd, strand, repName, repClass, repFamily \
                = (line.decode('utf-8').rstrip().split('\t')[i] for i in keep_fields)

            r_item = genome_generic(genoName, genoStart, genoEnd, strand)
            #will need to assign unique names to transcript_id
            dup_count = 1
            dup_num = transcript_id_unique.get(repName, 0)
            if dup_num < 1:
                unq_repName = repName 
                transcript_id_unique[repName] = 1
            else:
                unq_repName = repName + '_dup' + str(dup_num)
                transcript_id_unique[repName] += 1
                
                            
            r_item.add_meta({'score':score, 'gene_id':repName, 'transcript_id':unq_repName, 
                            'family_id':repFamily, 'class_id':repClass, 
                            'all_id':','.join([repName,repClass,repFamily,unq_repName])})
            # all_id = Repeat name (name, repName), 
            # Class (repeat_type_name, repClass), 
            # Family (repeat_subtype_name, repFamily), 
            # Other information (transcript name, accession)

            # print(r_item.to_gtf('mm10_rmsk', 'exon'))
            out_gtf.write(r_item.to_gtf('mm10_rmsk', 'exon', original_data_0_based=True))
            #     count += 1
            # else:
            #     break





if __name__ == "__main__":
    main()