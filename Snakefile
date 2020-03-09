'''
Author: Peter Chovanec
Aim: A Snakemake workflow to process RNA-DNA and DNA-DNA SPRITE-seq data
'''


import os 
import sys


#Location of scripts
barcode_id_jar = "scripts/java/BarcodeIdentification_v1.2.0.jar"
lig_eff = "scripts/python/get_ligation_efficiency.py"
split_fq = "scripts/split_dpm_rpm_fq.py"
atttb = "scripts/add_tnx_tag_to_bam.py"
add_chr = "scripts/python/ensembl2ucsc.py"
get_clusters = "scripts/python/get_clusters.py"
comb_anno = "scripts/combine_annotation_bams.py"
hicorrector = "scripts/HiCorrector_1.2/bin/ic"
clusters_heatmap = "scripts/python/get_sprite_contacts.py"

#Load config.yaml file

try:
    config_path = config["config_path"]
except:
    config_path = 'config.yaml'

configfile: config_path

try:
    email = config['email']
except:
    print("Won't send email on error")
    email = None

try:
    out_dir = config['output_dir']
    print('All data will be written to:', out_dir)
except:
    out_dir = ''
    print('Defaulting to working directory as output directory')

try:
    bid_config = config['bID']
    print('Using BarcodeID config', bid_config)
except:
    bid_config = 'workup/config.txt'
    print('Config "bID" not specified, looking for config at:', bid_config)

try:
    num_tags = config['num_tags']
    print('Using', num_tags, 'tags')
except:
    num_tags = "5"
    print('Config "num_tags" not specified, using:', num_tags)

#Make pipeline compatible for multiple assemblies
try:
    assembly = config['assembly']
    assert assembly in ['mm10', 'hg19', 'hg38'], 'Only "mm10" or "hg19" or "hg38" currently supported'
    print('Using', assembly)
except:
    print('Config "assembly" not specified, defaulting to "mm10"')
    assembly = 'mm10'

try:
    sprite_type = config['type']
    assert sprite_type in ['DNA-DNA', 'RNA-DNA'], 'Only "DNA-DNA" or "RNA-DNA" currently supported'
    print(sprite_type, 'SPRITE')
except:
    print('Config "type" not specified, defaulting to "RNA-DNA"')
    sprite_type = 'RNA-DNA'

try:
    anno_gtf = config['anno_gtf'][config['assembly']]
    anno_repeats_gtf = config['anno_repeats_gtf'][config['assembly']]
    mask = config['mask'][config['assembly']]
    exon_intron_gtf = config['exon_intron_gtf'][config['assembly']]
except:
    print('Annotation or mask path not specified in config.yaml')
    sys.exit() #no default, exit

try:
    run_snpsplit = config['snpsplit']
    print('Running SNPsplit:', run_snpsplit)
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
    print('SNPsplit not specified in config, will not run')
    snp_file = ''
    run_snpsplit = False
    snpsplit_name = ''
    anno_out_dir = 'alignments'


if run_snpsplit == 'True':
    try:
        g1 = config['snp_alleles'][config['cell']]['g1']
        print('Genome 1 strain:', g1)
    except:
        print('No strain specified for genome 1')
        g1 = None

    try:
        g2 = config['snp_alleles'][config['cell']]['g2']
        print('Genome 2 stain:', g2)
    except:
        print('No strain specified for genome 2')
        g2 = None    
else:
    g1 = None
    g2 = None    

try:
    DNA_aligner = config['dna_aligner']
    print('Using',DNA_aligner, 'for DNA alignment')
except:
    print('DNA aligner not specified, defaulting to bowtie2')
    DNA_aligner = 'bowtie2'

if run_snpsplit == 'True' and DNA_aligner == 'star':
        print('Running SNPsplit, will use Bowtie2')
        DNA_aligner = "bowtie2"


try:
    RNA_aligner = config['rna_aligner']
    if RNA_aligner == 'hisat2':
        hisat2_index = config['hisat2_index'][config['assembly']]
        hisat2_ss = config['hisat2_splice_sites'][config['assembly']]
    print('Using', RNA_aligner, 'for RNA alignment')
except:
    print('RNA aligner not specified in config.yaml')
    sys.exit()

try:
    if DNA_aligner == 'star' or RNA_aligner == 'star':
        star_index = config['star_index'][config['assembly']]
        star_repeat_index = config['star_repeat_index'][config['assembly']]
except:
    print('STAR indexes path not specified in config.yaml')
    sys.exit() #no default, exit

try:
    if DNA_aligner == 'bowtie2':
        if run_snpsplit == 'True':
            bowtie2_index = config['bowtie2_index'][config['cell']]
        else:
            bowtie2_index = config['bowtie2_index'][config['assembly']]
        bowtie2_repeat_index = config['bowtie2_repeat_index'][config['assembly']]
except:
    print('Bowtie2 index not specified in config.yaml')
    sys.exit() #no default, exit



#human reference has been generated with STAR 2.5.3a!

RNA_star_params = "--runMode alignReads \
--outFilterMultimapNmax 50 \
--outFilterScoreMinOverLread 0.30 \
--outFilterMatchNminOverLread 0.30 \
--outFilterIntronMotifs None \
--alignIntronMax 50000 \
--alignMatesGapMax 1000 \
--genomeLoad NoSharedMemory \
--outReadsUnmapped Fastx \
--alignIntronMin 80 \
--alignSJDBoverhangMin 5 \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 25000000000 \
--outSAMattributes All \
--readFilesCommand zcat \
--sjdbOverhang 100 \
--twopassMode Basic"
#--quantMode GeneCounts \

DNA_star_params = "--runMode alignReads \
--outFilterMultimapNmax 50 \
--outFilterScoreMinOverLread 0.30 \
--outFilterMatchNminOverLread 0.30 \
--outFilterIntronMotifs None \
--alignIntronMax 1 \
--alignMatesGapMax 2000 \
--genomeLoad NoSharedMemory \
--outReadsUnmapped Fastx \
--alignSJDBoverhangMin 5 \
--outSAMtype BAM Unsorted \
--limitBAMsortRAM 25000000000 \
--outSAMattributes All \
--readFilesCommand zcat \
--twopassMode Basic \
--alignEndsType EndToEnd"
#required by SNPsplit
# --outSAMtype BAM Unsorted
# --outSAMattributes NH HI NM MD
# --alignEndsType EndToEnd

#--sjdbOverhang 100 \
#--quantMode GeneCounts \

try:
    samples = config['samples']
    # print('Using samples file:', samples)
except:
    samples = './samples.json'
    # print('Defaulting to working directory for samples json file')

#get all samples from fastq Directory using the fastq2json.py scripts, then just
#load the json file with the samples
FILES = json.load(open(samples))
ALL_SAMPLES = sorted(FILES.keys())

ALL_FASTQ = []
for SAMPLE, file in FILES.items():
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R1')])
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R2')])

#Shared
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

#RNA-DNA
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

#STAR alignment
STAR_DNA_ALIGN = expand(out_dir + "workup/alignments/{sample}.DNA.Aligned.out.mapq20.bam",
                    sample=ALL_SAMPLES)
STAR_ALL_RNAr = expand(out_dir + "workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.bam", 
                       sample=ALL_SAMPLES)
STAR_RNA_UNMAP = expand(out_dir + "workup/alignments/{sample}.RNAr.unmapped.fastq.gz", 
                        sample=ALL_SAMPLES)
STAR_RNA = expand(out_dir + "workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.mapq20.bam", 
                   sample=ALL_SAMPLES)
STAR_RNAr = expand(out_dir + "workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.mapq20.bam", 
                  sample=ALL_SAMPLES)
STAR_TAG_ALL = expand(out_dir + "workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.mapq20.tag.bam", 
                 sample=ALL_SAMPLES)
STAR_ANNO_RNA = expand(out_dir + "workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.mapq20.bam.featureCounts.bam", 
                  sample=ALL_SAMPLES)



#DNA-DNA
BARCODEID_DNA = expand(out_dir + "workup/fastqs/{sample}_{read}.barcoded_dpm.fastq.gz", 
                       sample = ALL_SAMPLES, read = ["R1", "R2"])
BCS_DNA = expand(out_dir + "workup/alignments/{sample}.DNAonly.chr.bam", sample=ALL_SAMPLES)
CLUSTERS_DNA = expand(out_dir + "workup/clusters/{sample}.DNA.clusters", sample=ALL_SAMPLES)
MAPQ_DNA = expand(out_dir + "workup/alignments/{sample}.DNAonly.Aligned.out.mapq20.bam", sample=ALL_SAMPLES)




# if assembly == 'mm10':
if sprite_type == 'RNA-DNA':
    if DNA_aligner == 'star' and RNA_aligner == 'star':
        rule all:
            input: ALL_FASTQ + TRIM + TRIM_LOG + TRIM_RD + BARCODEID + LE_LOG_ALL + SPLIT_ALL +
                STAR_DNA_ALIGN + STAR_RNA + STAR_ALL_RNAr + STAR_ANNO_RNA +
                STAR_RNA_UNMAP + STAR_TAG_ALL + CHR_ALL + MASKED + 
                CLUSTERS + MULTI_QC
    elif DNA_aligner == 'bowtie2' and RNA_aligner == 'hisat2':
        rule all:
            input: ALL_FASTQ + TRIM + TRIM_LOG + TRIM_RD + BARCODEID + LE_LOG_ALL + SPLIT_ALL +
                Bt2_DNA_ALIGN + SNPSPLIT_DNA + Ht2_RNA_ALIGN + Bt2_TAG_ALL +  
                Ht2_ANNO_RNA + RNA_COMBINE + CHR_ALL + MASKED + CLUSTERS + MULTI_QC
    else:
        print('Combination not configured yet')
        sys.exit()
elif sprite_type == 'DNA-DNA':
    #just use bowtie2 here
    rule all:
        input: ALL_FASTQ + TRIM + TRIM_LOG + BARCODEID_DNA + LE_LOG_ALL + BARCODEID_DNA +
               Bt2_DNA_ALIGN + CHR_DNA + CLUSTERS_DNA + MULTI_QC



#Send and email if an error occurs during execution
onerror:
    shell('mail -s "an error occurred" ' + email + ' < {log}')


####################################################################################################
#Trimming and barcode identification
####################################################################################################

#Trim adaptors
#multiple cores requires pigz to be installed on the system
rule adaptor_trimming_pe:
    input:
        [lambda wildcards: FILES[wildcards.sample]['R1'],
        lambda wildcards: FILES[wildcards.sample]['R2']]
    output:
         out_dir + "workup/trimmed/{sample}_R1_val_1.fq.gz",
         out_dir + "workup/trimmed/{sample}_R1.fastq.gz_trimming_report.txt",
         out_dir + "workup/trimmed/{sample}_R2_val_2.fq.gz",
         out_dir + "workup/trimmed/{sample}_R2.fastq.gz_trimming_report.txt"
    log:
        out_dir + "workup/logs/{sample}.trim_galore.logs"
    conda:
        "envs/trim_galore.yaml"
    shell:
        "trim_galore \
        --paired \
        --gzip \
        --quality 20 \
        --fastqc \
        -o {out_dir}workup/trimmed/ \
        {input} &> {log}"



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
        adapters_r1 = "-a GATCGGAAGAG -a ATCAGCACTTA -g GGTGGTCTTT -g GCCTCTTGTT \
        -g CCAGGTATTT -g TAAGAGAGTT -g TTCTCCTCTT -g ACCCTCGATT",
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
        r1_barcoded = out_dir + "workup/fastqs/{sample}_R1.barcoded.fastq.gz" if sprite_type == 'RNA-DNA' 
        else out_dir + "workup/fastqs/{sample}_R1.barcoded_dpm.fastq.gz",
        r2_barcoded = out_dir + "workup/fastqs/{sample}_R2.barcoded.fastq.gz" if sprite_type == 'RNA-DNA'
        else out_dir + "workup/fastqs/{sample}_R2.barcoded_dpm.fastq.gz"
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


#TODO: this breaks the DNA-DNA pipeline
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
        "envs/bowtie2.yaml"
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
    

rule star_align_dna:
     input:
         fq = out_dir + "workup/fastqs/{sample}_R1.barcoded_dpm.fastq.gz"
     output:
         out_dir + "workup/alignments/{sample}.DNA.Aligned.out.mapq20.bam",
         out_dir + "workup/alignments/{sample}.DNA.Aligned.out.bam",
         out_dir + "workup/alignments/{sample}.DNA.Log.final.out",
         out_dir + "workup/alignments/{sample}.DNA.Log.out",
         out_dir + "workup/alignments/{sample}.DNA.Log.progress.out",
         out_dir + "workup/alignments/{sample}.DNA.SJ.out.tab",
         out_dir + "workup/alignments/{sample}.DNA.unmapped.fastq.gz"
     log:
         out_dir + "workup/logs/{sample}.DNA.star.log"
     threads: 10
     conda:
        'envs/star.yaml' if assembly == 'mm10' or assembly == 'hg38' else 
        'envs/star_hg19.yaml'
     shell:
         '''
         STAR {DNA_star_params} \
         --runThreadN {threads} \
         --genomeDir {star_index} \
         --readFilesIn {input.fq} \
         --outFileNamePrefix {out_dir}workup/alignments/{wildcards.sample}.DNA. &> {log}

         mv {out_dir}workup/alignments/{wildcards.sample}.DNA.Unmapped.out.mate1 \
             {out_dir}workup/alignments/{wildcards.sample}.DNA.unmapped.fastq

         pigz {out_dir}workup/alignments/{wildcards.sample}.DNA.unmapped.fastq

         samtools view -bq 20 {out_dir}workup/alignments/{wildcards.sample}.DNA.Aligned.out.bam > \
             {out_dir}workup/alignments/{wildcards.sample}.DNA.Aligned.out.mapq20.bam
         '''


#OPTIONAL: split this into 3 separate rules to make it run faster
rule annotate_dna:
    '''
    Users can specify the ‘-M’ option to fully count every alignment 
    reported for a multi-mapping read (each alignment carries 1 count.
    
    Users can specify the‘-O’ option to fully count them for each overlapping 
    meta-feature/feature (each overlapping meta-feature/feature 
    receives a count of 1 from a read (snoRNA's in introns of genes)
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
        -R BAM -M -s 1 \
        -g gene_name -a {exon_intron_gtf} -o {output.counts_exon} \
        {params.ex_rename}

        mv {params.ex_rename} {params.in_rename}
        featureCounts -T {threads} -t intron \
        -R BAM -M -s 1 -O \
        -g gene_name -a {exon_intron_gtf} -o {output.counts_intron} \
        {params.in_rename}

        mv {params.in_rename} {params.r_rename}
        featureCounts -T {threads} -t exon \
        -R BAM -M -s 1 \
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
        "envs/python_dep.yaml"
    shell:
        "python {comb_anno} -i {input.bam_exon} {input.bam_intron} {input.bam_rrna} \
            -i2 {input.bam} \
            -o {output} &> {log}"






############################################################################################
#RNA alignment
############################################################################################



rule star_align_rna:
    '''
    Align RNA with STAR to the genome first, annotate repeats, 
    anything that does not align, realign with bowtie2 to our repeats reference
    '''
    input:
         #  fq = "workup/alignments/{sample}.RNAr.unmapped.fastq.gz"
        fq = out_dir + "workup/fastqs/{sample}_R1.barcoded_rpm.fastq.gz"
    output:
        out_dir + "workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.mapq20.bam",
        out_dir + "workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.bam",
        out_dir + "workup/alignments/{sample}.RNA.Log.final.out",
        out_dir + "workup/alignments/{sample}.RNA.Log.out",
        out_dir + "workup/alignments/{sample}.RNA.Log.progress.out",
        out_dir + "workup/alignments/{sample}.RNA.SJ.out.tab",
        out_dir + "workup/alignments/{sample}.RNA.unmapped.fastq.gz"
    log:
        out_dir + "workup/logs/{sample}.RNA.star.log"
    threads: 10
    conda:
        'envs/star.yaml' if assembly == 'mm10' or assembly == 'hg38' else 
        'envs/star_hg19.yaml'
    shell:
        '''
        STAR {RNA_star_params} \
        --runThreadN {threads} \
        --genomeDir {star_index} \
        --readFilesIn {input.fq} \
        --outFileNamePrefix {out_dir}workup/alignments/{wildcards.sample}.RNA. &> {log}

        mv {out_dir}workup/alignments/{wildcards.sample}.RNA.Unmapped.out.mate1 \
            {out_dir}workup/alignments/{wildcards.sample}.RNA.unmapped.fastq

        pigz {out_dir}workup/alignments/{wildcards.sample}.RNA.unmapped.fastq

        samtools view -bq 20 {out_dir}workup/alignments/{wildcards.sample}.RNA.Aligned.sortedByCoord.out.bam > \
            {out_dir}workup/alignments/{wildcards.sample}.RNA.Aligned.sortedByCoord.out.mapq20.bam
        '''



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
        "envs/hisat2.yaml"
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
        "envs/python_dep.yaml"
    shell:
        "python {comb_anno} -i {input.bam_exon} {input.bam_intron} {input.bam_rrna} \
            -i2 {input.bam} \
            -o {output} &> {log}"



rule annotate_rna_star:
    input:
        out_dir + "workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.mapq20.bam"
    threads: 10
    output:
        bam=out_dir + "workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.mapq20.bam.featureCounts.bam",
        counts=out_dir + "workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.mapq20.bam.featureCounts.txt"
    log:
        out_dir + "workup/logs/{sample}.anno.log"
    conda:
        "envs/annotate_rna.yaml"
    shell:
        "featureCounts \
        -T {threads} \
        -t gene \
        -R BAM \
        -M \
        -s 1 \
        -g gene_name \
        -a {anno_gtf} \
        -o {output.counts} \
        {input}"


############################################################################################
#Repeats alignment
############################################################################################
#Align RNA with Bowtie2 to repeats
#MapQ filter 20, -F 4 only mapped reads, -F 256 remove not primary alignment reads
#-F: Do not output alignments with any bits set in INT present in the FLAG field
#q 20 don't filter until we resolve multimapping issue
rule bowtie2_align_rna_repeats:
     input:
        fq = out_dir + "workup/alignments/{sample}.RNA.hisat2.unmapped.lowmq.fq.gz"
     output:
         out_dir + "workup/alignments/{sample}.RNAr.bowtie2.mapq20.bam"   
     log:
         out_dir + "workup/logs/{sample}.RNAr.bowtie2.log"
     threads: 10
     conda:
          "envs/bowtie2.yaml"
     shell:
         "(bowtie2 \
        -p 10 \
        -t \
        --phred33 \
        -x {bowtie2_repeat_index} \
        -U {input.fq} | \
        samtools view -bS -F 4 -F 256 - > {output}) &> {log}"


#Align RNA with STAR to repeats
rule star_align_rrna:
     input:
         fq = out_dir + "workup/alignments/{sample}.RNA.unmapped.fastq.gz"
     output:
         out_dir + "workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.mapq20.bam",
         out_dir + "workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.bam",        
         out_dir + "workup/alignments/{sample}.RNAr.Log.final.out",
         out_dir + "workup/alignments/{sample}.RNAr.Log.out",
         out_dir + "workup/alignments/{sample}.RNAr.Log.progress.out",
         out_dir + "workup/alignments/{sample}.RNAr.SJ.out.tab",
         out_dir + "workup/alignments/{sample}.RNAr.unmapped.fastq.gz"
     log:
         out_dir + "workup/logs/{sample}.RNAr.star.log"
     threads: 10
     conda:
         'envs/star.yaml'
     shell:
         '''
         STAR {RNA_star_params} \
         --runThreadN {threads} \
         --genomeDir {star_repeat_index} \
         --readFilesIn {input.fq} \
         --outFileNamePrefix {out_dir}workup/alignments/{wildcards.sample}.RNAr. &> {log}

         mv {out_dir}workup/alignments/{wildcards.sample}.RNAr.Unmapped.out.mate1 workup/alignments/{wildcards.sample}.RNAr.unmapped.fastq

         gzip {out_dir}workup/alignments/{wildcards.sample}.RNAr.unmapped.fastq

         samtools view -bq 20 {out_dir}workup/alignments/{wildcards.sample}.RNAr.Aligned.sortedByCoord.out.bam > \
             {out_dir}workup/alignments/{wildcards.sample}.RNAr.Aligned.sortedByCoord.out.mapq20.bam
         '''


#Add chromosome of repeats as tag the featureCounts uses
rule add_tags_star:
    input:
        out_dir + "workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.mapq20.bam"
    output:
        out_dir + "workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.mapq20.tag.bam"
    log:
        out_dir + "workup/logs/{sample}.RNAr.star.tag.log"
    conda:
        'envs/rna_repeats.yaml'
    shell:'''
        samtools index {input}

        python {atttb} -i {input} -o {output}
        '''

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
        'envs/rna_repeats.yaml'
    shell:'''
        samtools sort -T {wildcards.sample} -@ {threads} -o {output.sort} {input}
        samtools index {output.sort}

        python {atttb} -i {output.sort} -o {output.out}
        '''


#TODO: remove DNA part of pipeline
############################################################################################
#DNA-DNA only 
############################################################################################


rule add_chr_DNA_DNA:
    input:
        out_dir + "workup/alignments/{sample}.DNA.bowtie2.mapq20.bam"
    output:
        out_dir + "workup/alignments/{sample}.DNAonly.chr.bam"
    log:
        out_dir + "workup/logs/{sample}.DNA_chr.log"
    conda:
        "envs/python_dep.yaml"
    shell:
        '''
        python {add_chr} -i {input} -o {output} --assembly {assembly} &> {log}
        '''


rule make_clusters_DNA:
    input:
        out_dir + "workup/alignments/{sample}.DNA.chr.masked.bam",
    output:
        out_dir + "workup/clusters/{sample}.DNA.clusters"
    log:
        out_dir + "workup/clusters/{sample}.make_clusters.log"
    conda:
        "envs/python_dep.yaml"
    shell:
        "python {get_clusters} -i {input} -o {output} -n {num_tags} &> {log}"

############################################################################################
#Shared
############################################################################################

rule repeat_mask:
    input:
        out_dir + "workup/alignments/{sample}.DNA.chr.bam" if sprite_type == 'RNA-DNA' else
        out_dir + "workup/alignments/{sample}.DNAonly.chr.bam"
    output:
        out_dir + "workup/alignments/{sample}.DNA.chr.masked.bam"
    conda:
        "envs/bedtools.yaml"
    shell:
        '''
        bedtools intersect -v -a {input} -b {mask} > {output}
        '''

############################################################################################
#RNA-DNA only 
############################################################################################

rule add_chr_RNA_DNA:
    input:
        dpm=out_dir + "workup/alignments/{sample}.DNA.Aligned.out.mapq20.bam" 
            if DNA_aligner == 'star' else out_dir + f"workup/{anno_out_dir}/{{sample}}.DNA.bowtie2.mapq20.anno.bam",
        rpm=out_dir + "workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.mapq20.bam.featureCounts.bam"
            if RNA_aligner == 'star' else out_dir + "workup/alignments/{sample}.RNA.hisat2.mapq20.anno.bam",
        rpm_repeat=out_dir + "workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.mapq20.tag.bam"
                   if RNA_aligner == 'star' else out_dir + "workup/alignments/{sample}.RNAr.bowtie2.mapq20.tag.bam"
    output:
        dpm=out_dir + "workup/alignments/{sample}.DNA.chr.bam",
        rpm=out_dir + "workup/alignments/{sample}.RNA.chr.bam",
        rpm_repeat=out_dir + "workup/alignments/{sample}.RNAr.chr.bam"
    log:
        dpm=out_dir + "workup/logs/{sample}.DNA_bcs.log",
        rpm=out_dir + "workup/logs/{sample}.RNA_bcs.log",
        rpm_repeat=out_dir + "workup/logs/{sample}.RNAr_bcs.log"
    conda:
        "envs/python_dep.yaml"
    shell:
        '''
        python {add_chr} -i {input.dpm} -o {output.dpm} --assembly {assembly} &> {log.dpm}
        python {add_chr} -i {input.rpm} -o {output.rpm} --assembly {assembly} &> {log.rpm}
        python {add_chr} -i {input.rpm_repeat} -o {output.rpm_repeat} --assembly none &> {log.rpm_repeat}
        '''



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
        "envs/python_dep.yaml"
    shell:
        "python {get_clusters} \
        -i {input.dpm} {input.rpm} {input.rpm_repeat} \
        -o {output} \
        -n {num_tags} \
        -g1 {g1} \
        -g2 {g2} &> {log}"




rule multiqc:
    input:
        #needs to be the last file produced in the pipeline 
        expand(out_dir + "workup/clusters/{sample}.clusters", sample=ALL_SAMPLES) if sprite_type=="RNA-DNA" else
        expand(out_dir + "workup/clusters/{sample}.DNA.clusters", sample=ALL_SAMPLES)
    output:
        out_dir + "workup/qc/multiqc_report.html"
    log:
        out_dir + "workup/logs/multiqc.log"
    conda: 
        "envs/qc.yaml"
    shell: 
        "multiqc {out_dir}workup -o {out_dir}workup/qc"




# rule snpsplit_DNA:
#     input:
#         "workup/alignments/{sample}.DNA.chr.bam"
#     output:
#         "workup/alignments/{sample}.DNA.chr.allele_flagged.bam"
#     log:
#         "workup/logs/snpsplit.log"
#     conda:
#         "envs/snpsplit.yaml"
#     shell:
#         "SNPsplit --snp_file {snp_file} {input} &> {log}"
