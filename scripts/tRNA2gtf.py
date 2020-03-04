import re
import gzip


#%%
trna = '/mnt/data/genomes/GRCm38.p6/repeats/mm10-tRNAs.fa'
trna_gtf = '/mnt/data/genomes/GRCm38.p6/repeats/mm10-tRNAs.gtf.gz'

#>Mus_musculus_tRNA-Ala-AGC-2-2 (tRNAscan-SE ID: chr13.trna89) Ala (AGC) 72 bp Sc: 84.7 chr13:21234045-21234116 (+)
with open(trna) as in_fa, gzip.open(trna_gtf, 'wt') as out_gtf:
    for qname, seq in fasta_parse(in_fa):
        name, *_, score, location, strand = qname.split(' ')
        gname = name.split('_')[-1]
        chrom, start, end = re.split(':|-', location)
        seqname = chrom.strip('chr')
        source = 'gtrnadb'
        feature = 'exon'
        frame = '.'
        strand = strand.strip('()')
        attribute = 'all_id ' + '"' + gname + '"'
        all_fields = [seqname, source, feature, start, end, score, strand, frame, attribute]
        out_gtf.write('\t'.join(all_fields) + '\n')

#%%
def fasta_parse(fp):
    '''Parse fasta files

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
        elif line.startswith('@'):
            raise Exception('Input starts with @ suggesting this is a fastq no fasta')
        else:
            seq += rec        