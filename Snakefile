'''
Author: Peter Chovanec
Aim: A Snakemake workflow to process RNA-DNA SPRITE-seq data
'''


import os 
import sys
import datetime
from pathlib import Path
import yaml

#when making dag.pdf all print statements need to be commented out, otherwise it will cause an error!

################################################################################
#Location of scripts
################################################################################

barcode_id_jar = "scripts/java/BarcodeIdentification_v1.2.0.jar"
lig_eff = "scripts/python/get_ligation_efficiency.py"
split_fq = "scripts/python/split_dpm_rpm_fq.py"
atttb = "scripts/python/add_tnx_tag_to_bam.py"
add_chr = "scripts/python/ensembl2ucsc.py"
get_clusters = "scripts/python/get_clusters.py"
comb_anno = "scripts/python/combine_annotation_bams.py"
hicorrector = "scripts/HiCorrector_1.2/bin/ic"
clusters_heatmap = "scripts/python/get_sprite_contacts.py"

################################################################################
#Load config.yaml file and other general settings
################################################################################

#Copy config file into logs
v = datetime.datetime.now()
run_date = v.strftime('%Y.%m.%d.')

# Print messages before start
INFO = ''


try:
    email = config['email']
except:
    INFO += "Won't send email on error \n"
    email = None

try:
    out_dir = config['output_dir']
    INFO += 'All data will be written to: ' + out_dir + '\n'
except:
    out_dir = ''
    INFO += 'Defaulting to working directory as output directory \n'

try:
    bid_config = config['bID']
    INFO += 'Using BarcodeID config ' + bid_config + '\n'
except:
    bid_config = 'workup/config.txt'
    INFO += 'Config "bID" not specified, looking for config at: ' + bid_config + '\n'

try:
    num_tags = config['num_tags']
    INFO += 'Using ' + num_tags + ' tags' + '\n'
except:
    num_tags = "5"
    INFO += 'Config "num_tags" not specified, using: ' + num_tags + '\n'


try:
    assembly = config['assembly']
    assert assembly in ['mm10', 'hg38'], 'Only "mm10" or "hg38" currently supported'
    INFO += 'Using ' + assembly + '\n'
except:
    sys.exit('ERROR: Config "assembly" not specified')


try:
    samples = config['samples']
    INFO += 'Using samples file: ' + samples + '\n'
except:
    samples = './samples.json'
    INFO += 'Defaulting to working directory for samples json file\n' 


################################################################################
#Annotation
################################################################################
try:
    anno_gtf = config['anno_gtf'][config['assembly']]
    anno_repeats_gtf = config['anno_repeats_gtf'][config['assembly']]
    mask = config['mask'][config['assembly']]
    exon_intron_gtf = config['exon_intron_gtf'][config['assembly']]
except:
    sys.exit('ERROR: Annotation or mask path not specified in config.yaml') #no default, exit

################################################################################
#SNPsplit 
################################################################################

try:
    run_snpsplit = config['snpsplit']
    INFO += 'Running SNPsplit: ' + run_snpsplit + '\n'
    if run_snpsplit == 'True':
        snp_file = config['snp_file'][config['cell']]
        snpsplit_name = '.allele_flagged'
        anno_out_dir = 'SNPsplit'
    else:
        snp_file = ''
        run_snpsplit = False
        snpsplit_name = ''
        anno_out_dir = 'alignments'
except:
    INFO += 'SNPsplit not specified in config, will not run\n' 
    snp_file = ''
    run_snpsplit = False
    snpsplit_name = ''
    anno_out_dir = 'alignments'


if run_snpsplit == 'True':
    try:
        g1 = config['snp_alleles'][config['cell']]['g1']
        INFO += 'Genome 1 strain: ' + g1 + '\n'
    except:
        INFO += 'No strain specified for genome 1\n'
        g1 = None

    try:
        g2 = config['snp_alleles'][config['cell']]['g2']
        INFO += 'Genome 2 stain: ' + g2 + '\n'
    except:
        INFO += 'No strain specified for genome 2\n'
        g2 = None    
else:
    g1 = None
    g2 = None    

################################################################################
#Aligner settings
################################################################################

try:
    hisat2_index = config['hisat2_index'][config['assembly']]
    hisat2_ss = config['hisat2_splice_sites'][config['assembly']]
except:
    sys.exit('ERROR: HISAT2 indexes not specified in config.yaml')


try:
    if run_snpsplit == 'True':
        bowtie2_index = config['bowtie2_index'][config['cell']]
    else:
        bowtie2_index = config['bowtie2_index'][config['assembly']]
    bowtie2_repeat_index = config['bowtie2_repeat_index'][config['assembly']]
except:
    sys.exit('ERROR: Bowtie2 index not specified in config.yaml') #no default, exit

################################################################################
#make output directories (aren't created automatically on cluster)
################################################################################

Path(out_dir + "workup/logs/cluster").mkdir(parents=True, exist_ok=True)
out_created = os.path.exists(out_dir + "workup/logs/cluster")
INFO += 'Output logs path created: ' + str(out_created) + '\n'

################################################################################
#Setup out files
################################################################################


#get all samples from fastq Directory using the fastq2json.py scripts, then just
#load the json file with the samples
FILES = json.load(open(samples))
ALL_SAMPLES = sorted(FILES.keys())

ALL_FASTQ = []
for SAMPLE, file in FILES.items():
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R1')])
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R2')])

CONFIG = [out_dir + "workup/logs/config_" + run_date + "yaml"]


TRIM = expand(out_dir + "workup/trimmed/{sample}_{read}.fq.gz", sample = ALL_SAMPLES, 
              read = ["R1_val_1", "R2_val_2"])
TRIM_LOG = expand(out_dir + "workup/trimmed/{sample}_{read}.fastq.gz_trimming_report.txt", 
                  sample = ALL_SAMPLES, read = ["R1", "R2"])
TRIM_RD = expand([out_dir + "workup/trimmed/{sample}_R1_val_1_RDtrim.fq.gz", 
                  out_dir + "workup/trimmed/{sample}_R2_val_2_RDtrim.fq.gz"], 
                  sample = ALL_SAMPLES)
LE_LOG_ALL = [out_dir + "workup/ligation_efficiency.txt"]
MASKED = expand(out_dir + "workup/alignments/{sample}.DNA.chr.masked.bam", sample=ALL_SAMPLES)
MULTI_QC = [out_dir + "workup/qc/multiqc_report.html"]

BARCODEID = expand(out_dir + "workup/fastqs/{sample}_{read}.barcoded.fastq.gz", sample = ALL_SAMPLES, 
                   read = ["R1", "R2"])
SPLIT_ALL = expand([out_dir + "workup/fastqs/{sample}_R1.barcoded_rpm.fastq.gz", 
                    out_dir + "workup/fastqs/{sample}_R1.barcoded_dpm.fastq.gz"], 
                    sample=ALL_SAMPLES)

CHR_ALL = expand([out_dir + "workup/alignments/{sample}.DNA.chr.bam",
          out_dir + "workup/alignments/{sample}.RNA.chr.bam",
          out_dir + "workup/alignments/{sample}.RNAr.chr.bam"], sample=ALL_SAMPLES)
RNA_COMBINE = expand(out_dir + "workup/alignments/{sample}.RNA.hisat2.mapq20.anno.bam", 
                     sample=ALL_SAMPLES)
CLUSTERS = expand(out_dir + "workup/clusters/{sample}.clusters", sample=ALL_SAMPLES)

#If aligning to N-masked genome for SNPsplit
#Bowtie2 alignment
Bt2_DNA_ALIGN = expand(out_dir + "workup/alignments/{sample}.DNA.bowtie2.mapq20.bam", 
                       sample=ALL_SAMPLES)
SNPSPLIT_DNA = expand(out_dir + f"workup/{anno_out_dir}/{{sample}}.DNA.bowtie2.mapq20{snpsplit_name}.bam", 
                      sample=ALL_SAMPLES)
Bt2_TAG_ALL = expand(out_dir + "workup/alignments/{sample}.RNAr.bowtie2.mapq20.tag.bam", 
                       sample=ALL_SAMPLES)
Bt2_RNAr = expand(out_dir + "workup/alignments/{sample}.RNAr.bowtie2.mapq20.bam", 
                  sample=ALL_SAMPLES)
DNA_COMBINE = expand(out_dir + f"workup/{anno_out_dir}/{{sample}}.DNA.bowtie2.mapq20.anno.bam", 
                     sample=ALL_SAMPLES)

#Hisat2 alignment
Ht2_RNA_ALIGN = expand([out_dir + "workup/alignments/{sample}.RNA.hisat2.mapq20.bam",
                        out_dir + "workup/alignments/{sample}.RNA.hisat2.unmapped.lowmq.fq.gz"], 
                        sample=ALL_SAMPLES)
Ht2_ANNO_RNA = expand([out_dir + "workup/alignments/{sample}.RNAex.hisat2.mapq20.bam.featureCounts.bam",
                       out_dir + "workup/alignments/{sample}.RNAin.hisat2.mapq20.bam.featureCounts.bam",
                      out_dir + "workup/alignments/{sample}.RNAr.hisat2.mapq20.bam.featureCounts.bam"],
                  sample=ALL_SAMPLES)


################################################################################
################################################################################
# Execute before workflow starts
################################################################################
################################################################################
onstart:
    print(INFO)

################################################################################
################################################################################
#Rule all
################################################################################
################################################################################

rule all:
    input: CONFIG + ALL_FASTQ + TRIM + TRIM_LOG + TRIM_RD + BARCODEID + LE_LOG_ALL + SPLIT_ALL +
        Bt2_DNA_ALIGN + SNPSPLIT_DNA + Ht2_RNA_ALIGN + Bt2_TAG_ALL +  
        Ht2_ANNO_RNA + RNA_COMBINE + CHR_ALL + MASKED + CLUSTERS + MULTI_QC



#Send and email if an error occurs during execution
onerror:
    shell('mail -s "an error occurred" ' + email + ' < {log}')


################################################################################
# Log
################################################################################

rule log_config:
    '''Copy config.yaml and place in logs folder with the date run
    '''
    output:
        out_dir + "workup/logs/config_" + run_date + "yaml"
    run:
        with open(output[0], 'w') as out:
            yaml.dump(config, out, default_flow_style=False)

################################################################################
#Merge
################################################################################
#For files that don't need to be merged, just create a symbolic link
rule merge_fastqs_pe:
    input: 
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output: 
        r1 = out_dir + "workup/merged/{sample}_R1.fastq.gz",
        r2 = out_dir + "workup/merged/{sample}_R2.fastq.gz"
    shell:
        ''' 
        count_1=$(echo '{input.r1}' | awk -F' ' '{{print NF}}')
        
        if [[ $count_1 -gt 1 ]]
        then
            cat {input.r1} > {output.r1}
        else
            ln -s {input.r1} {output.r1}
        fi

        count_2=$(echo '{input.r2}' | awk -F' ' '{{print NF}}')
        
        if [[ $count_2 -gt 1 ]]
        then
            cat {input.r2} > {output.r2}
        else
            ln -s {input.r2} {output.r2}
        fi
        '''

################################################################################
#Trimming and barcode identification
################################################################################

#Trim adaptors
#multiple cores requires pigz to be installed on the system
rule adaptor_trimming_pe:
    input:
        [out_dir + "workup/merged/{sample}_R1.fastq.gz", 
        out_dir + "workup/merged/{sample}_R2.fastq.gz"] 
        # [lambda wildcards: FILES[wildcards.sample]['R1'],
        # lambda wildcards: FILES[wildcards.sample]['R2']]
    output:
         out_dir + "workup/trimmed/{sample}_R1_val_1.fq.gz",
         out_dir + "workup/trimmed/{sample}_R1.fastq.gz_trimming_report.txt",
         out_dir + "workup/trimmed/{sample}_R2_val_2.fq.gz",
         out_dir + "workup/trimmed/{sample}_R2.fastq.gz_trimming_report.txt"
    threads:
        10
    log:
        out_dir + "workup/logs/{sample}.trim_galore.logs"
    conda:
        "envs/trim_galore.yaml"
    shell:
        '''
        if [[ {threads} -gt 8 ]]
        then 
            cores=2
        else
            cores=1
        fi

        trim_galore \
        --paired \
        --gzip \
        --cores $cores \
        --quality 20 \
        --fastqc \
        -o {out_dir}workup/trimmed/ \
        {input} &> {log}
        '''


rule cutadapt:
    '''
    Trim DPM RPM if read through reads
    TODO: Would be nice to run fastqc after this
    RPM from right ATCAGCACTTA
    DPM from right GATCGGAAGAG
    DPM from left GGTGGTCTT ^ anchored (only appears at the start of read)
    DPM5bot2-B1    /5Phos/TGACTTGTCATGTCTTCCGATCTGGTGGTCTTT
    DPM5bot3-C1    /5Phos/TGACTTGTCATGTCTTCCGATCTGCCTCTTGTT
    DPM5bot26-B4   /5Phos/TGACTTGTCATGTCTTCCGATCTCCAGGTATTT
    DPM5bot44-D6   /5Phos/TGACTTGTCATGTCTTCCGATCTTAAGAGAGTT
    DPM5bot85-E11  /5Phos/TGACTTGTCATGTCTTCCGATCTTTCTCCTCTT
    DPM5bot95-G12  /5Phos/TGACTTGTCATGTCTTCCGATCTACCCTCGATT
    '''
    input:
        [out_dir + "workup/trimmed/{sample}_R1_val_1.fq.gz", 
        out_dir + "workup/trimmed/{sample}_R2_val_2.fq.gz"]
    output:
        fastq1=out_dir + "workup/trimmed/{sample}_R1_val_1_RDtrim.fq.gz",
        fastq2=out_dir + "workup/trimmed/{sample}_R2_val_2_RDtrim.fq.gz",
        qc=out_dir + "workup/trimmed/{sample}.RDtrim.qc.txt"
    threads: 10
    params:
        adapters_r1 = "-a GATCGGAAGAG -a ATCAGCACTTA -g file:dpm96.fasta",
        adapters_r2 = "",
        others = "--minimum-length 20"
    log:
        "logs/cutadapt/{sample}.log"
    wrapper:
        "0.38.0/bio/cutadapt/pe"



#Identify barcodes using BarcodeIdentification_v1.2.0.jar
rule barcode_id:
    input:
        r1 = out_dir + "workup/trimmed/{sample}_R1_val_1_RDtrim.fq.gz",
        r2 = out_dir + "workup/trimmed/{sample}_R2_val_2_RDtrim.fq.gz"
    output:
    #if statements have to be inline (each input is like a function)
        r1_barcoded = out_dir + "workup/fastqs/{sample}_R1.barcoded.fastq.gz",
        r2_barcoded = out_dir + "workup/fastqs/{sample}_R2.barcoded.fastq.gz"
    log:
        out_dir + "workup/logs/{sample}.bID.log"
    shell:
        "java -jar {barcode_id_jar} \
        --input1 {input.r1} --input2 {input.r2} \
        --output1 {output.r1_barcoded} --output2 {output.r2_barcoded} \
        --config {bid_config} &> {log}"


#Get ligation efficiency
rule get_ligation_efficiency:
    input:
        r1 = out_dir + "workup/fastqs/{sample}_R1.barcoded.fastq.gz"
    output:
        temp(out_dir + "workup/{sample}.ligation_efficiency.txt")
    shell:
        "python {lig_eff} {input.r1} > {output}"


#Combine ligation efficiency from all samples into a single file
rule cat_ligation_efficiency:
    input:
        expand(out_dir + "workup/{sample}.ligation_efficiency.txt", sample=ALL_SAMPLES)
    output:
        out_dir + "workup/ligation_efficiency.txt"
    shell:
        "tail -n +1 {input} > {output}"



rule split_rpm_dpm:
    '''
    split rpm and dpm will also remove incomplete barcodes
    '''
    input:
        out_dir + "workup/fastqs/{sample}_R1.barcoded.fastq.gz"
    output:
        out_dir + "workup/fastqs/{sample}_R1.barcoded_dpm.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R1.barcoded_rpm.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R1.barcoded_other.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R1.barcoded_short.fastq.gz"
    log:
        out_dir + "workup/logs/{sample}_RPM_DPM.log"
    shell:
        "python {split_fq} --r1 {input} &> {log}"

############################################################################################
#DNA alignment
############################################################################################


rule bowtie2_align:
    '''
    MapQ filter 20, -F 4 only mapped reads, -F 256 remove not primary alignment reads
    -F: Do not output alignments with any bits set in INT present in the FLAG field
    '''
    input:
        fq=out_dir + "workup/fastqs/{sample}_R1.barcoded_dpm.fastq.gz"
    output:
        out_dir + "workup/alignments/{sample}.DNA.bowtie2.mapq20.bam"
    threads: 10
    log:
        out_dir + "workup/logs/{sample}.bowtie2.log"
    conda:
        "envs/alignment.yaml"
    shell:
        "(bowtie2 \
        -p 10 \
        -t \
        --phred33 \
        -x {bowtie2_index} \
        -U {input.fq} | \
        samtools view -bq 20 -F 4 -F 256 - > {output}) &> {log}"


rule snpsplit:
    input:
        out_dir + "workup/alignments/{sample}.DNA.bowtie2.mapq20.bam"
    output:
        out_dir + f"workup/{anno_out_dir}/{{sample}}.DNA.bowtie2.mapq20{snpsplit_name}.bam",
        out_dir + f"workup/{anno_out_dir}/{{sample}}.DNA.bowtie2.mapq20.genome1.bam",
        out_dir + f"workup/{anno_out_dir}/{{sample}}.DNA.bowtie2.mapq20.genome2.bam",
        out_dir + f"workup/{anno_out_dir}/{{sample}}.DNA.bowtie2.mapq20.SNPsplit_report.txt",
        out_dir + f"workup/{anno_out_dir}/{{sample}}.DNA.bowtie2.mapq20.SNPsplit_sort.txt",
        out_dir + f"workup/{anno_out_dir}/{{sample}}.DNA.bowtie2.mapq20.unassigned.bam"
    log:
        out_dir + "workup/logs/{sample}.snpsplit.log"
    conda:
        "envs/snpsplit.yaml"
    shell:
        "SNPsplit --snp_file {snp_file} -o {out_dir}workup/{anno_out_dir}/ {input} &> {log}"
    


rule annotate_dna:
    '''
    Users can specify the ‘-M’ option to fully count every alignment 
    reported for a multi-mapping read (each alignment carries 1 count.
    
    Users can specify the‘-O’ option to fully count them for each overlapping 
    meta-feature/feature (each overlapping meta-feature/feature 
    receives a count of 1 from a read (snoRNA's in introns of genes)

    Perform unstranded read counting on DNA
    '''
    input:
        out_dir + f"workup/{anno_out_dir}/{{sample}}.DNA.bowtie2.mapq20{snpsplit_name}.bam"
    threads: 10
    output:
        bam_exon=out_dir + f"workup/{anno_out_dir}/{{sample}}.DNAex.bowtie2.mapq20{snpsplit_name}.bam.featureCounts.bam",
        bam_intron=out_dir + f"workup/{anno_out_dir}/{{sample}}.DNAin.bowtie2.mapq20{snpsplit_name}.bam.featureCounts.bam",
        counts_exon=out_dir + f"workup/{anno_out_dir}/{{sample}}.DNAex.bowtie2.mapq20{snpsplit_name}.bam.featureCounts.txt",
        counts_intron=out_dir + f"workup/{anno_out_dir}/{{sample}}.DNAin.bowtie2.mapq20{snpsplit_name}.bam.featureCounts.txt",
        bam_rrna=out_dir + f"workup/{anno_out_dir}/{{sample}}.DNAr.bowtie2.mapq20{snpsplit_name}.bam.featureCounts.bam",
        counts_rrna=out_dir + f"workup/{anno_out_dir}/{{sample}}.DNAr.bowtie2.mapq20{snpsplit_name}.bam.featureCounts.txt"
    log:
        out_dir + "workup/logs/{sample}.anno.log"
    params:
        ex_rename = out_dir + f"workup/{anno_out_dir}/{{sample}}.DNAex.bowtie2.mapq20{snpsplit_name}.bam",
        in_rename = out_dir + f"workup/{anno_out_dir}/{{sample}}.DNAin.bowtie2.mapq20{snpsplit_name}.bam",
        r_rename = out_dir + f"workup/{anno_out_dir}/{{sample}}.DNAr.bowtie2.mapq20{snpsplit_name}.bam"
    conda:
        "envs/annotate_rna.yaml"
    shell:
        '''
        mv {input} {params.ex_rename}
        featureCounts -T {threads} -t exon \
        -R BAM -M -s 0 \
        -g gene_name -a {exon_intron_gtf} -o {output.counts_exon} \
        {params.ex_rename}

        mv {params.ex_rename} {params.in_rename}
        featureCounts -T {threads} -t intron \
        -R BAM -M -s 0 -O \
        -g gene_name -a {exon_intron_gtf} -o {output.counts_intron} \
        {params.in_rename}

        mv {params.in_rename} {params.r_rename}
        featureCounts -T {threads} -t exon \
        -R BAM -M -s 0 \
        -g all_id -a {anno_repeats_gtf} -o {output.counts_rrna} \
        {params.r_rename}

        mv {params.r_rename} {input}

        '''

rule combine_annotations_dna:
    input:
        bam=out_dir + f"workup/{anno_out_dir}/{{sample}}.DNA.bowtie2.mapq20{snpsplit_name}.bam",
        bam_exon=out_dir + f"workup/{anno_out_dir}/{{sample}}.DNAex.bowtie2.mapq20{snpsplit_name}.bam.featureCounts.bam",
        bam_intron=out_dir + f"workup/{anno_out_dir}/{{sample}}.DNAin.bowtie2.mapq20{snpsplit_name}.bam.featureCounts.bam",
        bam_rrna=out_dir + f"workup/{anno_out_dir}/{{sample}}.DNAr.bowtie2.mapq20{snpsplit_name}.bam.featureCounts.bam"
    output:
        out_dir + f"workup/{anno_out_dir}/{{sample}}.DNA.bowtie2.mapq20.anno.bam"
    log:
        out_dir + "workup/logs/{sample}.DNA_anno_combine.log"
    conda:
        "envs/alignment.yaml"
    shell:
        "python {comb_anno} -i {input.bam_exon} {input.bam_intron} {input.bam_rrna} \
            -i2 {input.bam} \
            -o {output} &> {log}"




############################################################################################
#RNA alignment
############################################################################################

rule hisat2_align:
    '''
    #from clusterflow pipeline
    # we are currently using a very high penalty score for soft-clipping (--sp 1000,1000)
    #because Hisat2 seems to soft-clip even when it should run in --end-to-end mode
    # we are also filtering out unmapped reads (-F 4), or reads where the mate was unmapped (-F 8)
    # we are also filtering non-primary alignments (-F 256)
    #filter on mapq score of 20 (Skip alignments with MAPQ smaller than 20)
    #-U FILE Write alignments that are not selected by the various filter options to FILE
    '''
    input:
        fq=out_dir + "workup/fastqs/{sample}_R1.barcoded_rpm.fastq.gz"
    output:
        all_reads=temp(out_dir + "workup/alignments/{sample}.RNA.hisat2.bam"),
        # low_mapq=temp("workup/alignments/{sample}.RNA.hisat2.lowmapq.bam"),
        # unmapped=temp("workup/alignments/{sample}.RNA.hisat2.unmapped.bam"),
        mapped=out_dir + "workup/alignments/{sample}.RNA.hisat2.mapq20.bam",
        merged=out_dir + "workup/alignments/{sample}.RNA.hisat2.unmapped.lowmq.bam",
        fq_gz=out_dir + "workup/alignments/{sample}.RNA.hisat2.unmapped.lowmq.fq.gz"
    threads: 10
    conda:
        "envs/alignment.yaml"
    log:
        out_dir + "workup/logs/{sample}.hisat2.log"
    shell:
        '''
        (hisat2 --sp 1000,1000 \
        -p 10 \
        -t \
        --phred33 \
        --known-splicesite-infile {hisat2_ss} \
        -x {hisat2_index} \
        -U {input.fq} | \
        samtools view -b -F 256 - > {output.all_reads}) &> {log}
        #split out unmapped and low mapq reads for realignment to repeats
        samtools view -bq 20 -U {output.merged} -F 4 {output.all_reads} > {output.mapped}
        samtools bam2fq -@ {threads} {output.merged} > {out_dir}workup/alignments/{wildcards.sample}.RNA.hisat2.unmapped.lowmq.fq
        pigz {out_dir}workup/alignments/{wildcards.sample}.RNA.hisat2.unmapped.lowmq.fq
        '''


################################################################################
#RNA annotation
################################################################################


rule annotate_rna:
    '''
    -M                  Multi-mapping reads will also be counted. For a multi-
                        mapping read, all its reported alignments will be 
                        counted. The 'NH' tag in BAM/SAM input is used to detect 
                        multi-mapping reads.
    -s <int or string>  Perform strand-specific read counting. A single integer
                        value (applied to all input files) or a string of comma-
                        separated values (applied to each corresponding input
                        file) should be provided. Possible values include:
                        0 (unstranded), 1 (stranded) and 2 (reversely stranded).
                        Default value is 0 (ie. unstranded read counting carried
                        out for all input files).
    -t <string>         Specify feature type in GTF annotation. 'exon' by 
                        default. Features used for read counting will be 
                        extracted from annotation using the provided value.
    -g <string>         Specify attribute type in GTF annotation. 'gene_id' by 
                        default. Meta-features used for read counting will be 
                        extracted from annotation using the provided value.
        '''
    input:
        out_dir + "workup/alignments/{sample}.RNA.hisat2.mapq20.bam"
    threads: 10
    output:
        bam_exon=out_dir + "workup/alignments/{sample}.RNAex.hisat2.mapq20.bam.featureCounts.bam",
        bam_intron=out_dir + "workup/alignments/{sample}.RNAin.hisat2.mapq20.bam.featureCounts.bam",
        counts_exon=out_dir + "workup/alignments/{sample}.RNAex.hisat2.mapq20.bam.featureCounts.txt",
        counts_intron=out_dir + "workup/alignments/{sample}.RNAin.hisat2.mapq20.bam.featureCounts.txt",
        bam_rrna=out_dir + "workup/alignments/{sample}.RNAr.hisat2.mapq20.bam.featureCounts.bam",
        counts_rrna=out_dir + "workup/alignments/{sample}.RNAr.hisat2.mapq20.bam.featureCounts.txt"
    log:
        out_dir + "workup/logs/{sample}.anno.log"
    params:
        ex_rename = out_dir + "workup/alignments/{sample}.RNAex.hisat2.mapq20.bam",
        in_rename = out_dir + "workup/alignments/{sample}.RNAin.hisat2.mapq20.bam",
        r_rename = out_dir + "workup/alignments/{sample}.RNAr.hisat2.mapq20.bam"
    conda:
        "envs/annotate_rna.yaml"
    shell:
        '''
        mv {input} {params.ex_rename}
        featureCounts -T {threads} -t exon \
        -R BAM -M -s 1 \
        -g gene_name -a {exon_intron_gtf} -o {output.counts_exon} \
        {params.ex_rename}

        mv {params.ex_rename} {params.in_rename}
        featureCounts -T {threads} -t intron \
        -R BAM -M -s 1 \
        -g gene_name -a {exon_intron_gtf} -o {output.counts_intron} \
        {params.in_rename}

        mv {params.in_rename} {params.r_rename}
        featureCounts -T {threads} -t exon \
        -R BAM -M -s 1 \
        -g all_id -a {anno_repeats_gtf} -o {output.counts_rrna} \
        {params.r_rename}

        mv {params.r_rename} {input}

        '''


rule combine_annotations_rna:
    input:
        bam=out_dir + "workup/alignments/{sample}.RNA.hisat2.mapq20.bam",
        bam_exon=out_dir + "workup/alignments/{sample}.RNAex.hisat2.mapq20.bam.featureCounts.bam",
        bam_intron=out_dir + "workup/alignments/{sample}.RNAin.hisat2.mapq20.bam.featureCounts.bam",
        bam_rrna=out_dir + "workup/alignments/{sample}.RNAr.hisat2.mapq20.bam.featureCounts.bam"
    output:
        out_dir + "workup/alignments/{sample}.RNA.hisat2.mapq20.anno.bam"
    log:
        out_dir + "workup/logs/{sample}.RNA_anno_combine.log"
    conda:
        "envs/alignment.yaml"
    shell:
        "python {comb_anno} -i {input.bam_exon} {input.bam_intron} {input.bam_rrna} \
            -i2 {input.bam} \
            -o {output} &> {log}"





############################################################################################
#Repeats alignment
############################################################################################

rule bowtie2_align_rna_repeats:
    '''
    Align RNA with Bowtie2 to repeats
    MapQ filter 20, -F 4 only mapped reads, -F 256 remove not primary alignment reads
    -F: Do not output alignments with any bits set in INT present in the FLAG field
    q 20 don't filter until we resolve multimapping issue
    '''
    input:
        fq = out_dir + "workup/alignments/{sample}.RNA.hisat2.unmapped.lowmq.fq.gz"
    output:
        out_dir + "workup/alignments/{sample}.RNAr.bowtie2.mapq20.bam"   
    log:
        out_dir + "workup/logs/{sample}.RNAr.bowtie2.log"
    threads: 10
    conda:
        "envs/alignment.yaml"
    shell:
        "(bowtie2 \
        -p 10 \
        -t \
        --phred33 \
        -x {bowtie2_repeat_index} \
        -U {input.fq} | \
        samtools view -bS -F 4 -F 256 - > {output}) &> {log}"

#Add chromosome of repeats as tag the featureCounts uses
rule add_tags_bowtie2:
    input:
        out_dir + "workup/alignments/{sample}.RNAr.bowtie2.mapq20.bam"
    output:
        out=out_dir + "workup/alignments/{sample}.RNAr.bowtie2.mapq20.tag.bam",
        idx=temp(out_dir + "workup/alignments/{sample}.RNAr.bowtie2.mapq20.sorted.bam.bai"),
        sort=temp(out_dir + "workup/alignments/{sample}.RNAr.bowtie2.mapq20.sorted.bam")
    log:
        out_dir + "workup/logs/{sample}.RNAr.bowtie2.tag.log"
    threads:
        5
    conda:
        'envs/alignment.yaml'
    shell:'''
        samtools sort -T {wildcards.sample} -@ {threads} -o {output.sort} {input}
        samtools index {output.sort}

        python {atttb} -i {output.sort} -o {output.out}
        '''


############################################################################################
#Mask and filters
############################################################################################

rule repeat_mask:
    input:
        out_dir + "workup/alignments/{sample}.DNA.chr.bam"
    output:
        out_dir + "workup/alignments/{sample}.DNA.chr.masked.bam"
    conda:
        "envs/bedtools.yaml"
    shell:
        '''
        bedtools intersect -v -a {input} -b {mask} > {output}
        '''

rule add_chr:
    input:
        dpm=out_dir + f"workup/{anno_out_dir}/{{sample}}.DNA.bowtie2.mapq20.anno.bam",
        rpm=out_dir + "workup/alignments/{sample}.RNA.hisat2.mapq20.anno.bam",
        rpm_repeat=out_dir + "workup/alignments/{sample}.RNAr.bowtie2.mapq20.tag.bam"
    output:
        dpm=out_dir + "workup/alignments/{sample}.DNA.chr.bam",
        rpm=out_dir + "workup/alignments/{sample}.RNA.chr.bam",
        rpm_repeat=out_dir + "workup/alignments/{sample}.RNAr.chr.bam"
    log:
        dpm=out_dir + "workup/logs/{sample}.DNA_bcs.log",
        rpm=out_dir + "workup/logs/{sample}.RNA_bcs.log",
        rpm_repeat=out_dir + "workup/logs/{sample}.RNAr_bcs.log"
    conda:
        "envs/alignment.yaml"
    shell:
        '''
        python {add_chr} -i {input.dpm} -o {output.dpm} --assembly {assembly} &> {log.dpm}
        python {add_chr} -i {input.rpm} -o {output.rpm} --assembly {assembly} &> {log.rpm}
        python {add_chr} -i {input.rpm_repeat} -o {output.rpm_repeat} --assembly none &> {log.rpm_repeat}
        '''


################################################################################
#Make clusters
################################################################################


rule make_merged_clusters:
    input:
        dpm=out_dir + "workup/alignments/{sample}.DNA.chr.masked.bam",
        rpm=out_dir + "workup/alignments/{sample}.RNA.chr.bam",
        rpm_repeat=out_dir + "workup/alignments/{sample}.RNAr.chr.bam"
    output:
        out_dir + "workup/clusters/{sample}.clusters"
    log:
        out_dir + "workup/clusters/{sample}.make_clusters.log"
    conda:
        "envs/alignment.yaml"
    shell:
        "python {get_clusters} \
        -i {input.dpm} {input.rpm} {input.rpm_repeat} \
        -o {output} \
        -n {num_tags} \
        -g1 {g1} \
        -g2 {g2} &> {log}"


################################################################################
#MultiQC
################################################################################

rule multiqc:
    input:
        #needs to be the last file produced in the pipeline 
        expand(out_dir + "workup/clusters/{sample}.clusters", sample=ALL_SAMPLES)
    output:
        out_dir + "workup/qc/multiqc_report.html"
    log:
        out_dir + "workup/logs/multiqc.log"
    conda: 
        "envs/qc.yaml"
    shell: 
        "multiqc {out_dir}workup -o {out_dir}workup/qc"

