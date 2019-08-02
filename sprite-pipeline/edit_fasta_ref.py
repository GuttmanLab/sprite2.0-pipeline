
# from babrahamlinkon.general import fasta_parse

ncrna = '/mnt/data/genomes/custom/ncRNAs/ribornaandsmallrnas.fa'

new_ncrna = '/mnt/data/genomes/custom/ncRNAs/ribornaandsmallrnas_simple.fa'

total_len = 0
with open(ncrna,'r') as fa, \
open(new_ncrna, 'w') as out_fa:
    for qname, seq in fasta_parse(fa):
        # try:
        sp_name = qname.strip('>').split('|')
        # except:
            # sp_name = qname.strip('>')

        try:
            gene_name = sp_name[-2]
            # gene_biotype = sp_name[-1]
        except IndexError:
            gene_name = sp_name[0]
            # gene_biotype = ''

        #write out new interleaved fasta (80 characters per line)

        out_fa.write('>' + gene_name + '\n')
        len_out = 0
        initial_len = len(seq)
        total_len += initial_len
        while len_out < initial_len:
            out_fa.write(seq[:80] + '\n')
            seq = seq[80:]
            len_out += 80
print(total_len) #use for STAR genome size


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
