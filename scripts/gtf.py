from collections import defaultdict
import gzip
import re
from tqdm import tqdm

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


def meta2gtf(meta_data):
    '''Format a dictionary into gtf metadata format
    
    Args:
        meta_data(dict): dictionary with metadata
    '''
    out_f = []
    for k, v in meta_data.items():
        out_f.append(str(k) + ' ' + '"' + str(v) + '"')
    out = ';'.join(out_f)
    return out


def gtf2meta(meta_data):
    '''Format gtf metadata into a dictionary

    Args:
        meta_data(str): gtf metadata string
    '''
    meta_dict = defaultdict()
    pattern = re.compile('([a-zA-Z0-9_]+)\\ \"([a-zA-Z0-9\\ ]+)')
    fields = meta_data.split('; ')
    for i in fields:
        match = pattern.search(i)
        try:
            name, value = match.groups()
        except AttributeError:
            print(fields)
        meta_dict[name] = value

    return(meta_dict)


class genome_generic:
    '''A generic class to hold genomic information

    Methods:
    - add_meta(): Add a dictionary with additional metadata
    "{nameOfFeature:'myFeature'}"
    - to_gtf(): Return a GTF formated string
    '''

    def __init__(self, chrom, start, end, strand, f_type=None):
        self._chrom = chrom
        self._start = start
        self._end = end
        self._stand = strand
        self._type = f_type
        self._meta = defaultdict(set)

    def add_meta(self, metadata):
        if isinstance(metadata, dict):
            for k, v in metadata.items():
                #each feature name can only have a single value
                self._meta[k] = v

    def to_gtf(self, name, gtype, original_data_0_based=True):
        '''
        convert to 1-based if original data is 0-based
        '''
        if original_data_0_based:
            start = int(self._start) + 1
        else:
            start = int(self._start)
        
        out_format = '\t'.join([self._chrom, name, gtype, str(start), self._end, 
        self._meta.get('score', '0'), self._stand, '.', meta2gtf(self._meta), '\n'])
        return out_format

    def fetch_meta(self, field):
        '''
        fetch a specific metadata field
        '''
        out = self._meta.get(field, None)
        if out == None:
            avail_fields = ' '.join(list(self._meta.keys()))
            assert out != None, f'Field not present in metadata, present field: {avail_fields}'
        else:
            return out



def parse_gtf(path, feature_type=None):
    '''
    parse gtf file

    GFF is 1-based coordinate system
    gtf format:
    1       mm10_rmsk       exon    3044870 3045464 4839    +       .       gene_id "RLTR1D"; transcript_id "RLTR1D"; family_id "ERV1"; class_id "LTR";
    1       mm10_rmsk       exon    3014472 3014749 1604    +       .       gene_id "L1VL1"; transcript_id "L1VL1"; family_id "L1"; class_id "LINE"; all_id "L1VL1,L1VL1,L1,LINE"

    Args:
        path(str): Path of GTF file
        feature_type(str): Only parse: gene, transcript, exon, CDS, start_codon,
                           stop_codon, five_prime_utr, three_prime_utr, or 
                           Selenocysteine

    Return:
        list of genome_generic objects
    '''
    lines = list()
    with file_open(path) as gtf:
        for line in tqdm(gtf):
            line = line.decode('utf-8').rstrip().strip('\n')
            if not line.startswith('#'):
                chrom, _, f_type, start, end, _, strand, _, meta = line.split('\t')
                if feature_type != None:
                    if f_type == feature_type:
                        record = genome_generic(chrom, start, end, strand, f_type)
                        record.add_meta(gtf2meta(meta))
                        lines.append(record)
                else:
                    record = genome_generic(chrom, start, end, strand)
                    record.add_meta(gtf2meta(meta))
                    lines.append(record)

    return lines
