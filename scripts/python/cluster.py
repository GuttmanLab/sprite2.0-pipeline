import pysam
import re
import gzip
import os
import sys
from collections import defaultdict, Counter
from tqdm import tqdm

class Position:
    """This class represents a genomic position, with type of nucleic acid (RNA or DNA)

    Methods:
    - to_string(): Returns a string representation of this position in the form
      "R/DPM(feature)_chrX:1000"
    """

    def __init__(self, read_type, strand, feature, chromosome, start_coordinate, end_coordinate):
        self._type = read_type
        self._strand = strand
        self._feature = feature
        self._chromosome = chromosome
        self._start_coordinate = start_coordinate
        self._end_coordinate = end_coordinate

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False

        return (self._type == other._type and
                self._strand == other._strand and
                self._feature == other._feature and
                self._chromosome == other._chromosome and
                self._start_coordinate == other._start_coordinate and
                self._end_coordinate == other._end_coordinate)

    def __hash__(self):
        return hash((self._type, self._strand, self._feature, self._chromosome, 
                     self._start_coordinate, self._end_coordinate))

    def to_string(self):
        try:
            out = self._type + "[" + self._strand + ';' + self._feature + "]" \
                + "_" + self._chromosome + ":" + \
                str(self._start_coordinate) + "-" + str(self._end_coordinate)
        except:
            print(self._type, self._feature, self._chromosome)
            print('Elements are not as expect!')
            sys.exit()
        return out

class UMIs:
    """This class hold the UMI

    Methods:
        to_string: Returns a string of the UMI
    """

    def __init__(self):
        self._umi = []

    def add_umi(self, umi):
        self._umi.append(umi)

    def unique(self):
        return '\t'.join(set(self._umi))

    def to_string(self):
        return '\t'.join(self._umi)

class Cluster:
    """This class represents a barcoding cluster as a collection of genomic
    positions.

    The underlying data structure is a set, so duplicate positions are
    discarded.

    Methods:
    - add_position(position): Adds a genomic position to this cluster

    - size(): Returns the number of reads or positions in this cluster

    - to_string(): Returns a string representation of this cluster as a
      tab-delimtited series of positions. See Position#to_string for how
      positions are represented as strings.

    - to_list(): Returns the Position class as a list (similar to_string())
    """

    def __init__(self):
        self._positions = set()

    def __iter__(self):
        return iter(self._positions)

    def add_position(self, position):
        self._positions.add(position)

    def size(self, read_type=None):
        if read_type == None:
            return len(self._positions)
        else:
            return sum([1 if pos._type == read_type else 0 for pos in self._positions])

    def to_string(self):
        position_strings = [position.to_string() for position in self._positions]
        return "\t".join(position_strings)

    def to_list(self):
        position_strings = [position.to_string() for position in self._positions]
        return position_strings

    def count_type(self):
        rna_dna_count = Counter([position._type for position in self._positions])
        return rna_dna_count


class Clusters:
    """This class represents a collection of barcoding clusters.

    Methods:
    - get_cluster(barcode): Returns the cluster that corresponds to the given
      barcode. If the cluster does not exist, it is initialized (with zero
      positions), and this empty cluster is returned.

    - get_items(): iterate over clusters dictionary yielding keys and values

    - add_position(barcode, position): Adds the position to the cluster
      that corresponds with the given barcodes

    - to_strings(): Returns an iterator over the string representations of all
      of the contained clusters.

    - remove_cluster(barcode): Removes a cluster with the specified barcode

    - add_umi(barcode, umi): Used for fastq to cluster format, will used umi
      instead of position coordinates

    - unique(): keep only unique cluster entries

    - make_lookup(): make a lookup table for converting cluster back into bam
    """
    def __init__(self):
        self._clusters = {}

    def __iter__(self):
        return iter(self._clusters.values())

    def __getitem__(self, barcode):
        return self._clusters[barcode]

    def get_cluster(self, barcode):
        if barcode not in self._clusters:
            self._clusters[barcode] = Cluster()
        return self._clusters[barcode]

    def get_items(self):
        return self._clusters.items()
    
    def add_position(self, barcode, position):
        self.get_cluster(barcode).add_position(position)

    def to_strings(self):
        for barcode, cluster in self._clusters.items():
            yield barcode + '\t' + cluster.to_string()

    def remove_cluster(self, barcode):
        del self._clusters[barcode]

    def get_umi_cluster(self, barcode):
        if barcode not in self._clusters:
            self._clusters[barcode] = UMIs()
        return self._clusters[barcode]

    def add_umi(self, barcode, umi):
        self.get_umi_cluster(barcode).add_umi(umi)

    def unique(self):
        for barcode, cluster in self._clusters.items():
            yield barcode + '\t' + cluster.unique()

    def make_lookup(self):
        lookup = defaultdict(set)
        for barcode, cluster in self._clusters.items():
            lookup[barcode].update(cluster.to_list())
        return lookup


# test = Clusters()
# test.add_position('ewafa', Position('DPM', 'blabla', 'chr1', '123', '132',))
# for k, v in test.get_items():
#     print(k, v)
# for i in test:
#     print(i)
# test['ewafa']
# test['ewafa'].count_type()

def anno_type(element):
    '''
    Extract annotation type
    '''
    return element.split('.')[-1]


def order_annotation(anno):
    '''
    Enforce exon, intron, repeat order
    '''
    anno_list = anno.split(';')
    sorted_anno = sorted(anno_list, key=anno_type)
    return ';'.join(sorted_anno)


def get_clusters(bamfile, num_tags, genome_1, genome_2):
    """Parses a BAM file, groups positions into clusters according to their
    barcodes, and returns the resulting structure.

    Each BAM record must have the barcode stored in the query name like so:

    ORIGINAL_READ_NAME::[Tag1][Tag2][Tag3]

    The tags should be enclosed in brackets and separated from
    the original read name with a double-colon.
    """
    #strip RPM DPM from barcode
    #TODO add file name as an additional barcode
    
    clusters = Clusters()
    pattern = re.compile('::' + num_tags * '\[([a-zA-Z0-9_\-]+)\]')

    for bam in bamfile:
        #get sample name from bamfile
        file_name = os.path.basename(bam)
        sample_name = file_name.split('.')[0]
        try:
            with pysam.AlignmentFile(bam, "rb") as f:
                for read in f.fetch(until_eof = True):
                    name = read.query_name
                    match = pattern.search(name)
                    barcode = list(match.groups())
                    strand = '+' if not read.is_reverse else '-'
                    #strip RPM DPM from barcode
                    if 'RPM' in barcode:
                        #get featureCounts annotation
                        if read.has_tag('XT'):
                            gene_anno = read.get_tag('XT')
                        elif read.has_tag('XS'):
                            gene_anno = read.get_tag('XS')
                        else:
                            raise Exception('XS tag missing, was featureCounts run?')
                        
                        #[Strand; Gene annotation]
                        position = Position('RPM', strand, order_annotation(gene_anno), 
                                            read.reference_name, read.reference_start, 
                                            read.reference_end)
                        barcode.remove('RPM')
                    elif 'DPM' in barcode:
                        #Allele annotation
                        # XX:Z:G1
                        if read.has_tag('XX'):
                            allele = read.get_tag('XX')
                            if genome_1 != "None" and genome_2 != "None":
                                if allele == "G1":
                                    allele = genome_1
                                elif allele == "G2":
                                    allele = genome_2
                            allele_anno = allele + ';' + strand 
                        else:
                            allele_anno = ''
                        #Gene annotation
                        if read.has_tag('XT'):
                            gene_anno = read.get_tag('XT')
                        elif read.has_tag('XS'):
                            gene_anno = read.get_tag('XS')
                        else:
                            gene_anno = ''

                        #anno = Strand; Gene annotation; Allele
                        anno = order_annotation(gene_anno) + ';' + allele_anno

                        position = Position('DPM', strand, anno, read.reference_name,
                                            read.reference_start, read.reference_end)
                        barcode.remove('DPM')
                    else:
                        raise Exception('RPM or DPM not present in full barcode')
                    barcode.append(sample_name)
                    barcode_str = ".".join(barcode)
                    clusters.add_position(barcode_str, position)
        except ValueError:
            print('BAM file provided is not a BAM or is empty!')

    return clusters


def get_clusters_fastq(fastqfile, num_tags, orientation, umi_len):
    """Parses a fastq file, groups positions into clusters according to their
    barcodes, and returns the resulting structure.

    Each fastq record must have the barcode stored in the query name like so:

    ORIGINAL_READ_NAME::[Tag1][Tag2][Tag3]

    The tags should be enclosed in brackets and separated from
    the original read name with a double-colon.
    """
    #strip RPM DPM from barcode
    clusters = Clusters()
    pattern = re.compile('::' + num_tags * '\[([a-zA-Z0-9_\-]+)\]')
    count = 0
    with file_open(fastqfile) as f:
        for name, seq, thrd, qual in fastq_parse(f):
            #extract UMI and place as position data
            match = pattern.search(name)
            barcode = list(match.groups())
            if orientation == 'r1': #read 1 orientation, UMI on left
                umi = seq[:umi_len]
                # umi = seq
            elif orientation == 'r2': #read 2 orientation, UMI on right
                umi = seq[-umi_len:]
                # umi = seq
            else:
                raise Exception('Incorrect orientation specified')

            barcode_str = "|".join(barcode[:num_tags])
            clusters.add_umi(barcode_str, umi)

            count += 1
    print('Reads parsed:', count)
    return clusters

# umi_len = 16
# num_tags = 8
# fastqfile = '/mnt/data/20190704_10_1000_spritezero/assembled/complex_10_S2_L001_R1_001_val_1.assembled_rv_bID.fq.gz'
# test = get_clusters_fastq(fastqfile, num_tags, 'r2', umi_len)
#
# test.get_umi_cluster('NOT_FOUND|R5-2_bot_A2|NOT_FOUND|NOT_FOUND|NOT_FOUND|NOT_FOUND|NOT_FOUND|NOT_FOUND').to_string()



def write_clusters_to_file(clusters, outfile, unique=False):
    """Writes a Clusters object to a file"""

    count = 0
    with open(outfile, 'w') as f:

        if unique:
            for cluster_string in clusters.unique():
                f.write(cluster_string)
                f.write("\n")
                count += 1
        else:
            for cluster_string in clusters.to_strings():
                f.write(cluster_string)
                f.write("\n")
                count += 1
    print('Number of clusters written:',count)


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



def parse_cluster(c_file):
    '''
    Parse cluster file

    Args:
        c_file(str): input path of cluster file
    '''

    total_reads = 0
    clusters = Clusters()
    # pattern = re.compile('([a-zA-Z0-9]+)\[([a-zA-Z0-9_;\:\,\-\+\.\(\)\?]+)\]_([a-zA-Z0-9_\:\-\.\,\(\)]+):([0-9]+)\-([0-9]+)')
    pattern = re.compile('([a-zA-Z0-9]+)\[(.+)\]_(.+):([0-9]+)\-([0-9]+)')
    # match = pattern.search('DPM[UA;-;Ralbp1.intron;Unassigned_NoFeatures.none]_chr17:65883728-65883812')
    
    with file_open(c_file) as c:
        for line in tqdm(c):

            barcode, *reads = line.decode('utf-8').rstrip('\n').split('\t')
            for read in reads:
                total_reads += 1
                try:
                    match = pattern.search(read)
                    read_type, feature, chrom, start, end = match.groups()
                    #get strand from annotation
                    strand, *anno = feature.split(';')
                    anno = order_annotation(';'.join(anno))
                    
                    if strand == '+' or strand == '-':
                        position = Position(read_type, strand, anno, chrom, 
                                            start, end)
                    else:
                        position = Position(read_type, '.', order_annotation(feature), 
                                            chrom, start, end)
                    clusters.add_position(barcode, position)
                except:
                    print(read)
                    raise Exception('Pattern did not match above printed string')
    print('Total cluster reads:', total_reads)
    return(clusters)




def write_bam(cluster, num_tags, original_bam, output_bam, genome_1, genome_2):
    '''From a cluster make/subset bam file
    If barcode, chrom, start and end coordinates match, write out read into new BAM

    Args:
        cluster(Clusters): 
    '''

    #get sample name from bamfile
    file_name = os.path.basename(original_bam)
    sample_name = file_name.split('.')[0]

    pattern = re.compile('::' + num_tags * '\[([a-zA-Z0-9_\-]+)\]')

    #get genome build
    if 'hg38' in file_name:
        assembly = 'hg38'
    elif 'mm10' in file_name:
        assembly = 'mm10'

    read_lookup = cluster.make_lookup()
    out_reads = 0
    with pysam.AlignmentFile(original_bam, "rb") as in_bam:
        out_bam = pysam.AlignmentFile(output_bam, "wb", template = in_bam)
        for read in in_bam.fetch(until_eof = True):
            name = read.query_name
            match = pattern.search(name)
            barcode = list(match.groups())
            strand = '+' if not read.is_reverse else '-'
            if 'RPM' in barcode:
                #get featureCounts annotation
                if read.has_tag('XT'):
                    gene_anno = read.get_tag('XT')
                elif read.has_tag('XS'):
                    gene_anno = read.get_tag('XS')
                else:
                    raise Exception('XS tag missing, was featureCounts run?')

                position = Position('RPM', strand, order_annotation(gene_anno), read.reference_name,
                                    read.reference_start, read.reference_end)
                barcode.remove('RPM')

            elif 'DPM' in barcode:
                #Allele annotation
                # XX:Z:G1
                if read.has_tag('XX'):
                    allele = read.get_tag('XX')
                    if genome_1 != "None" and genome_2 != "None":
                        if allele == "G1":
                            allele = genome_1
                        elif allele == "G2":
                            allele = genome_2
                    allele_anno = allele + ';' + strand 
                else:
                    allele_anno = ''
                #Gene annotation
                if read.has_tag('XT'):
                    gene_anno = read.get_tag('XT')
                elif read.has_tag('XS'):
                    gene_anno = read.get_tag('XS')
                else:
                    gene_anno = ''

                #anno = Strand; Gene annotation; Allele
                anno = order_annotation(gene_anno) + ';' + allele_anno

                position = Position('DPM', strand, anno, read.reference_name,
                                    read.reference_start, read.reference_end)
                barcode.remove('DPM')
            else:
                raise Exception('RPM or DPM not present in full barcode')  
        
            barcode.append(sample_name)
            barcode_str = ".".join(barcode)

            barcode_reads = read_lookup.get(barcode_str, "Not present")
         
            if position.to_string() in barcode_reads:
                out_bam.write(read)
                out_reads += 1 

    out_bam.close()
    print('Total reads written:', out_reads)




