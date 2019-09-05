# samples = ["2p5DNA1","5pMerge1","5pRNA1"]
email = "pchovanec@lncrna.caltech.edu"


#Location of scripts
barcode_id_jar = "sprite-pipeline/java/BarcodeIdentification_v1.2.0.jar"
lig_eff = "sprite-pipeline/python/get_ligation_efficiency.py"
split_fq = "sprite-pipeline/split_dpm_rpm_fq.py"
atttb = "sprite-pipeline/add_tnx_tag_to_bam.py"
all_tags = "sprite-pipeline/python/filter_all_tags.py"
get_clusters = "sprite-pipeline/python/get_clusters.py"

#Load config.yaml file

try:
    config_path = config["config_path"]
except:
    config_path = "config.yaml"

configfile: config_path


try:
    bid_config = config['bID']
    print('Using BarcodeID config', bid_config)
except:
    bid_config = "workup/config.txt"
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
    assert assembly in ['mm10', 'hg19'], 'Only "mm10" or "hg19" currently supported'
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
    star_index = config['star_index'][config['assembly']]
    star_repeat_index = config['star_repeat_index'][config['assembly']]
    anno_gtf = config['anno_gtf'][config['assembly']]
    mask = config['mask'][config['assembly']]
except:
    print('STAR indexes or annotation, mask paths not specified in config.yaml')

try:
    snp_file = config['snp_file'][config['assembly']]


# if assembly == 'mm10':
#     star_index = "/groups/guttman/Peter/genomes/GRCm38.p6/star"
#     star_repeat_index = "/groups/guttman/Peter/genomes/ncRNA/star"
#     hisat2_index = "/groups/guttman/Peter/genomes/GRCm38.p6"
#     bowtie2_index = "/groups/guttman/Peter/genomes/GRCm38.p6"
#     anno_gtf = "/groups/guttman/Peter/genomes/GRCm38.p6/GRCm38.p6.annotation.gtf.gz"
#     mask = "/groups/guttman/Peter/genomes/GRCm38.p6/blacklist.rmsk.mm10.milliDivLessThan140.bed"
# elif assembly == 'hg19':
#     star_index = "/groups/guttman/genomes/homo_sapiens/hg19/STAR/hg19"
#     star_repeat_index = "/groups/guttman/dmariani/Hs-45S"
#     anno_gtf = "/groups/guttman/genomes/homo_sapiens/hg19/STAR/hg19/gencode.v19.annotation.gtf"
#     mask = "/groups/guttman/genomes/combinations/hg19_mm9/masks/hg19-mm9.gatk35-and-rmsk140.bed"

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




# rscript = "/central/software/R/3.5.0/bin/Rscript"
# cluster_plot = "/groups/guttman/software/sprite-pipeline/r/get_cluster_size_distribution.r"


#import json
#import snakemake
import os


#get all samples from fastq Directory using the fastq2json.py scripts, then just
#load the json file with the samples
FILES = json.load(open("./samples.json"))
ALL_SAMPLES = sorted(FILES.keys())

ALL_FASTQ = []
for SAMPLE, file in FILES.items():
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R1')])
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R2')])

# ALL_FASTQ = expand("raw/{sample}_{read}.fastq.gz", sample = ALL_SAMPLES, read = ["R1", "R2"])

#Shared
TRIM = expand("workup/trimmed/{sample}_{read}.fq.gz", sample = ALL_SAMPLES, read = ["R1_val_1", "R2_val_2"])
TRIM_LOG = expand("workup/trimmed/{sample}_{read}.fastq.gz_trimming_report.txt", sample = ALL_SAMPLES, read = ["R1", "R2"])
LE_LOG_ALL = ["workup/ligation_efficiency.txt"]
STAR_ALL_DNA = expand("workup/alignments/{sample}.DNA.Aligned.sortedByCoord.out.bam", sample=ALL_SAMPLES)
MASKED = expand("workup/alignments/{sample}.DNA.all_bcs.masked.bam", sample=ALL_SAMPLES)
MULTI_QC = ["workup/qc/multiqc_report.html"]

#RNA-DNA
BARCODEID = expand("workup/fastqs/{sample}_{read}.barcoded.fastq.gz", sample = ALL_SAMPLES, read = ["R1", "R2"])
RPM_ALL = expand("workup/fastqs/{sample}_R1.barcoded_rpm.fastq.gz", sample=ALL_SAMPLES)
DPM_ALL = expand("workup/fastqs/{sample}_R1.barcoded_dpm.fastq.gz", sample=ALL_SAMPLES)
STAR_ALL_RNAr = expand("workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.bam", sample=ALL_SAMPLES)
STAR_RNA_UNMAP = expand("workup/alignments/{sample}.RNAr.unmapped.fastq.gz", sample=ALL_SAMPLES)
STAR_ALL_RNA = expand("workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.bam", sample=ALL_SAMPLES)
MAPQ_ALL = expand(["workup/alignments/{sample}.DNA.Aligned.sortedByCoord.out.mapq20.bam",
           "workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.mapq20.bam",
           "workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.mapq20.bam"], sample=ALL_SAMPLES)
BCS_ALL = expand(["workup/alignments/{sample}.DNA.all_bcs.bam",
          "workup/alignments/{sample}.RNA.all_bcs.bam",
          "workup/alignments/{sample}.RNAr.all_bcs.bam"], sample=ALL_SAMPLES)
TAG_ALL = expand("workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.mapq20.tag.bam", sample=ALL_SAMPLES)
ANNO_RNA = expand("workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.mapq20.bam.featureCounts.bam", sample=ALL_SAMPLES)
CLUSTERS = expand("workup/clusters/{sample}.clusters", sample=ALL_SAMPLES)

#DNA-DNA
BARCODEID_DNA = expand("workup/fastqs/{sample}_{read}.barcoded_dpm.fastq.gz", sample = ALL_SAMPLES, read = ["R1", "R2"])
MAPQ_DNA = expand("workup/alignments/{sample}.DNA.Aligned.sortedByCoord.out.mapq20.bam", sample=ALL_SAMPLES)
BCS_DNA = expand("workup/alignments/{sample}.DNAonly.all_bcs.bam", sample=ALL_SAMPLES)
CLUSTERS_DNA = expand("workup/clusters/{sample}.DNA.clusters", sample=ALL_SAMPLES)
MAPQ_DNA = expand("workup/alignments/{sample}.DNAonly.Aligned.sortedByCoord.out.mapq20.bam", sample=ALL_SAMPLES)



# HISAT2_ALL = expand("workup/alignments/{sample}.hisat2.bam", sample=ALL_SAMPLES)
# BOWTIE2_ALL = expand("workup/alignments/{sample}.bowtie2.bam", sample=ALL_SAMPLES)
# MERGED_ALL = expand("workup/alignments/{sample}.merged.bam", sample=ALL_SAMPLES)
# MERGED_ALL_STAR = expand("workup/alignments/{sample}.merged.bam", sample=ALL_SAMPLES)
# OTHER_ALL = expand("workup/fastqs/{sample}_R1.barcoded_dpm.fastq.gz", sample=ALL_SAMPLES)


if assembly == 'mm10':
    if sprite_type == 'RNA-DNA':
        rule all:
            input: ALL_FASTQ + TRIM + TRIM_LOG + BARCODEID + LE_LOG_ALL + RPM_ALL +
                DPM_ALL + STAR_ALL_RNAr + STAR_ALL_DNA + STAR_ALL_RNA + 
                STAR_RNA_UNMAP + MAPQ_ALL + TAG_ALL + BCS_ALL + ANNO_RNA + MASKED + 
                CLUSTERS + MULTI_QC #+ OTHER_ALL
    elif sprite_type == 'DNA-DNA':
        rule all:
            input: ALL_FASTQ + TRIM + TRIM_LOG + BARCODEID_DNA + LE_LOG_ALL + BARCODEID_DNA +
                   CLUSTERS_DNA + STAR_ALL_DNA + MULTI_QC + MAPQ_DNA
elif assembly == 'hg19':
    if sprite_type == 'RNA-DNA':
        rule all:
            input: ALL_FASTQ + TRIM + TRIM_LOG + BARCODEID + LE_LOG_ALL + RPM_ALL +
                DPM_ALL + STAR_ALL_RNAr + STAR_ALL_DNA + STAR_ALL_RNA + 
                STAR_RNA_UNMAP + MAPQ_ALL + TAG_ALL + BCS_ALL + ANNO_RNA + MASKED + 
                CLUSTERS + MULTI_QC
    elif sprite_type == 'DNA-DNA':
        rule all:
            input: ALL_FASTQ + TRIM + TRIM_LOG + BARCODEID_DNA + LE_LOG_ALL + BARCODEID_DNA +
                   CLUSTERS_DNA + STAR_ALL_DNA + MULTI_QC + MAPQ_DNA

#Send and email if an error occurs during execution
onerror:
    shell('mail -s "an error occurred" ' + email + ' < {log}')


####################################################################################################
#Rules
####################################################################################################

#Trim adaptors
#multiple cores requires pigz to be installed on the system
rule adaptor_trimming_pe:
    input:
        [lambda wildcards: FILES[wildcards.sample]['R1'],
        lambda wildcards: FILES[wildcards.sample]['R2']]
    output:
         "workup/trimmed/{sample}_R1_val_1.fq.gz",
         "workup/trimmed/{sample}_R1.fastq.gz_trimming_report.txt",
         "workup/trimmed/{sample}_R2_val_2.fq.gz",
         "workup/trimmed/{sample}_R2.fastq.gz_trimming_report.txt"
    log:
        "workup/logs/{sample}.trim_galore.logs"
    conda:
        "envs/trim_galore.yaml"
    shell:
        "trim_galore \
        --paired \
        --gzip \
        --quality 20 \
        --fastqc \
        -o workup/trimmed/ \
        {input} &> {log}"


#Identify barcodes using BarcodeIdentification_v1.2.0.jar
rule barcode_id:
    input:
        r1 = "workup/trimmed/{sample}_R1_val_1.fq.gz",
        r2 = "workup/trimmed/{sample}_R2_val_2.fq.gz",
    output:
    #if statements have to be inline (each input is like a function)
        r1_barcoded = "workup/fastqs/{sample}_R1.barcoded.fastq.gz" if sprite_type == 'RNA-DNA' 
        else "workup/fastqs/{sample}_R1.barcoded_dpm.fastq.gz",
        r2_barcoded = "workup/fastqs/{sample}_R2.barcoded.fastq.gz" if sprite_type == 'RNA-DNA'
        else "workup/fastqs/{sample}_R2.barcoded_dpm.fastq.gz"
    log:
        "workup/logs/{sample}.bID.log"
    shell:
        "java -jar {barcode_id_jar} \
        --input1 {input.r1} --input2 {input.r2} \
        --output1 {output.r1_barcoded} --output2 {output.r2_barcoded} \
        --config {bid_config} &> {log}"


#Get ligation efficiency
rule get_ligation_efficiency:
    input:
        r1 = "workup/fastqs/{sample}_R1.barcoded.fastq.gz"
    output:
        temp("workup/{sample}.ligation_efficiency.txt")
    shell:
        "python {lig_eff} {input.r1} > {output}"


#Combine ligation efficiency from all samples into a single file
rule cat_ligation_efficiency:
    input:
        expand("workup/{sample}.ligation_efficiency.txt", sample=ALL_SAMPLES)
    output:
        "workup/ligation_efficiency.txt"
    shell:
        "tail -n +1 {input} > {output}"



rule split_rpm_dpm:
    input:
        "workup/fastqs/{sample}_R1.barcoded.fastq.gz"
    output:
        "workup/fastqs/{sample}_R1.barcoded_dpm.fastq.gz",
        "workup/fastqs/{sample}_R1.barcoded_rpm.fastq.gz",
        "workup/fastqs/{sample}_R1.barcoded_other.fastq.gz"
    log:
        "workup/logs/{sample}_RPM_DPM.log"
    shell:
        "python {split_fq} --r1 {input} &> {log}"





#Align RNA with STAR to repeats
rule star_align_rna_repeats:
     input:
         fq = "workup/fastqs/{sample}_R1.barcoded_rpm.fastq.gz"
     output:
         "workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.bam",        
         "workup/alignments/{sample}.RNAr.Log.final.out",
         "workup/alignments/{sample}.RNAr.Log.out",
         "workup/alignments/{sample}.RNAr.Log.progress.out",
         "workup/alignments/{sample}.RNAr.SJ.out.tab",
         "workup/alignments/{sample}.RNAr.unmapped.fastq.gz"
     log:
         "workup/logs/{sample}.RNAr.star.log"
     threads: 10
     conda:
         'envs/rna_repeats.yaml'
     shell:
         '''
         STAR {RNA_star_params} \
         --runThreadN {threads} \
         --genomeDir {star_repeat_index} \
         --readFilesIn {input.fq} \
         --outFileNamePrefix workup/alignments/{wildcards.sample}.RNAr. &> {log}

         mv workup/alignments/{wildcards.sample}.RNAr.Unmapped.out.mate1 workup/alignments/{wildcards.sample}.RNAr.unmapped.fastq

         gzip workup/alignments/{wildcards.sample}.RNAr.unmapped.fastq

         '''

#Align RNA with STAR
rule star_align_rna:
     input:
         # fq = "workup/fastqs/{sample}_R1.barcoded_rpm.fastq.gz"
         fq = "workup/alignments/{sample}.RNAr.unmapped.fastq.gz"
     output:
         "workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.bam",
         "workup/alignments/{sample}.RNA.Log.final.out",
         "workup/alignments/{sample}.RNA.Log.out",
         "workup/alignments/{sample}.RNA.Log.progress.out",
         "workup/alignments/{sample}.RNA.SJ.out.tab",
         "workup/alignments/{sample}.RNA.unmapped.fastq.gz"
     log:
         "workup/logs/{sample}.RNA.star.log"
     threads: 10
     conda:
        'envs/star.yaml' if assembly == 'mm10' else 'envs/star_hg19.yaml'
     shell:
         '''
         STAR {RNA_star_params} \
         --runThreadN {threads} \
         --genomeDir {star_index} \
         --readFilesIn {input.fq} \
         --outFileNamePrefix workup/alignments/{wildcards.sample}.RNA. &> {log}

         mv workup/alignments/{wildcards.sample}.RNA.Unmapped.out.mate1 workup/alignments/{wildcards.sample}.RNA.unmapped.fastq

         gzip workup/alignments/{wildcards.sample}.RNA.unmapped.fastq
         '''





rule star_align_dna:
     input:
         fq = "workup/fastqs/{sample}_R1.barcoded_dpm.fastq.gz"
     output:
         "workup/alignments/{sample}.DNA.Aligned.sortedByCoord.out.bam",
         "workup/alignments/{sample}.DNA.Log.final.out",
         "workup/alignments/{sample}.DNA.Log.out",
         "workup/alignments/{sample}.DNA.Log.progress.out",
         "workup/alignments/{sample}.DNA.SJ.out.tab",
         "workup/alignments/{sample}.DNA.unmapped.fastq.gz"
     log:
         "workup/logs/{sample}.DNA.star.log"
     threads: 10
     conda:
        'envs/star.yaml'
     shell:
         """
         STAR {DNA_star_params} \
         --runThreadN {threads} \
         --genomeDir {star_index} \
         --readFilesIn {input.fq} \
         --outFileNamePrefix workup/alignments/{wildcards.sample}.DNA. &> {log}

         mv workup/alignments/{wildcards.sample}.DNA.Unmapped.out.mate1 workup/alignments/{wildcards.sample}.DNA.unmapped.fastq

         gzip workup/alignments/{wildcards.sample}.DNA.unmapped.fastq
         """



rule mapq_filter_RNA_DNA:
    input:
        dpm="workup/alignments/{sample}.DNA.Aligned.sortedByCoord.out.bam",
        rpm="workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.bam",
        rpm_repeat="workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.bam"
    output:
        dpm="workup/alignments/{sample}.DNA.Aligned.sortedByCoord.out.mapq20.bam",
        rpm="workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.mapq20.bam",
        rpm_repeat="workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.mapq20.bam"
    log:
        "workup/logs/{sample}.DNA_mapq.log"
        "workup/logs/{sample}.RNA_mapq.log"
        "workup/logs/{sample}.RNAr_mapq.log"
    conda:
        "envs/samtools.yaml"
    shell:
        '''
        samtools view -b -q 20 -o {output.dpm} {input.dpm}
        samtools view -b -q 20 -o {output.rpm} {input.rpm}
        samtools view -b -q 20 -o {output.rpm_repeat} {input.rpm_repeat}
        '''

rule mapq_filter_DNA_DNA:
    input:
        dpm="workup/alignments/{sample}.DNA.Aligned.sortedByCoord.out.bam"
    output:
        dpm="workup/alignments/{sample}.DNAonly.Aligned.sortedByCoord.out.mapq20.bam"
    log:
        "workup/logs/{sample}.DNA_mapq.log"
    conda:
        "envs/samtools.yaml"
    shell:
        '''
        samtools view -b -q 20 -o {output.dpm} {input.dpm}
        '''


rule add_tags:
    input:
        "workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.mapq20.bam"
    output:
        "workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.mapq20.tag.bam"
    log:
        "workup/logs/{sample}.RNAr.tag.log"
    conda:
        'envs/rna_repeats.yaml'
    shell:'''
        samtools index {input}

        python {atttb} -i {input} -o {output}
        '''


rule annotate_rna:
    input:
        "workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.mapq20.bam"
    threads: 8
    output:
        bam="workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.mapq20.bam.featureCounts.bam",
        counts="workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.mapq20.bam.featureCounts.txt"
        # "workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.featureCounts.bam.summary"
    log:
        "workup/logs/{sample}.anno.log"
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



#merge RNA DNA bam files
# rule merge_bam_star:
#     input:
#         dpm="workup/alignments/{sample}.DNA.Aligned.sortedByCoord.out.bam",
#         rpm="workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.bam.featureCounts.bam",
#         rpm_repeat="workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.tag.bam"
#     threads: 10
#     output:
#         dpm=temp("workup/alignments/{sample}.tempDPM.unique.bam"),
# 	    # rpm=temp("workup/alignments/{sample}.tempRPM.unique.bam"),
# 	    merged="workup/alignments/{sample}.merged.bam"
#     log:
#         "workup/logs/{sample}.merge.log"
#     conda:
#         "envs/samtools.yaml"
#     shell:'''
#         {samtools} view -b -q 255 -o {output.dpm} {input.dpm}
#         {samtools} merge -@ {threads} {output.merged} {output.dpm} {input.rpm} {input.rpm_repeat} &> {log}
#         '''

rule all_barcodes_RNA_DNA:
    input:
        # "workup/alignments/{sample}.merged.bam"
        dpm="workup/alignments/{sample}.DNA.Aligned.sortedByCoord.out.mapq20.bam",
        rpm="workup/alignments/{sample}.RNA.Aligned.sortedByCoord.out.mapq20.bam.featureCounts.bam",
        rpm_repeat="workup/alignments/{sample}.RNAr.Aligned.sortedByCoord.out.mapq20.tag.bam"
    output:
        # "workup/alignments/{sample}.merged.all_bcs.bam"
        dpm="workup/alignments/{sample}.DNA.all_bcs.bam",
        rpm="workup/alignments/{sample}.RNA.all_bcs.bam",
        rpm_repeat="workup/alignments/{sample}.RNAr.all_bcs.bam"
    log:
        # "workup/logs/{sample}.all_bcs.log"
        dpm="workup/logs/{sample}.DNA_bcs.log",
        rpm="workup/logs/{sample}.RNA_bcs.log",
        rpm_repeat="workup/logs/{sample}.RNAr_bcs.log"
    conda:
        "envs/python_dep.yaml"
    shell:
        '''
        python {all_tags} -i {input.dpm} -o {output.dpm} --assembly {assembly} &> {log.dpm}
        python {all_tags} -i {input.rpm} -o {output.rpm} --assembly {assembly} &> {log.rpm}
        python {all_tags} -i {input.rpm_repeat} -o {output.rpm_repeat} --assembly none &> {log.rpm_repeat}
        '''


rule all_barcodes_DNA_DNA:
    input:
        "workup/alignments/{sample}.DNAonly.Aligned.sortedByCoord.out.mapq20.bam",
    output:
        "workup/alignments/{sample}.DNAonly.all_bcs.bam",
    log:
        "workup/logs/{sample}.DNA_bcs.log",
    conda:
        "envs/python_dep.yaml"
    shell:
        '''
        python {all_tags} -i {input} -o {output} --assembly {assembly} &> {log}
        '''




rule repeat_mask:
    input:
        # "workup/alignments/{sample}.merged.all_bcs.bam"{samtools} view -b -q 255 -o {output.dpm} {input.dpm}
        "workup/alignments/{sample}.DNA.all_bcs.bam" if sprite_type == 'RNA-DPM' else
        "workup/alignments/{sample}.DNAonly.all_bcs.bam"
    output:
        "workup/alignments/{sample}.DNA.all_bcs.masked.bam"
    conda:
        "envs/bedtools.yaml"
    shell:
        '''
        bedtools intersect -v -a {input} -b {mask} > {output}
        '''

rule make_merged_clusters:
    input:
        dpm="workup/alignments/{sample}.DNA.all_bcs.masked.bam",
        rpm="workup/alignments/{sample}.RNA.all_bcs.bam",
        rpm_repeat="workup/alignments/{sample}.RNAr.all_bcs.bam"
    output:
        "workup/clusters/{sample}.clusters"
    log:
        "workup/clusters/{sample}.make_clusters.log"
    conda:
        "envs/python_dep.yaml"
    shell:
        "python {get_clusters} -i {input.dpm} {input.rpm} {input.rpm_repeat} -o {output} -n {num_tags} &> {log}"



rule make_clusters_DNA:
    input:
        "workup/alignments/{sample}.DNA.all_bcs.masked.bam",
    output:
        "workup/clusters/{sample}.DNA.clusters"
    log:
        "workup/clusters/{sample}.make_clusters.log"
    conda:
        "envs/python_dep.yaml"
    shell:
        "python {get_clusters} -i {input} -o {output} -n {num_tags} &> {log}"




rule multiqc:
    input:
        #needs to be the last file produced in the pipeline 
        expand("workup/clusters/{sample}.clusters", sample=ALL_SAMPLES) if sprite_type=="RNA-DNA" else
        expand("workup/clusters/{sample}.DNA.clusters", sample=ALL_SAMPLES)
    output:
        "workup/qc/multiqc_report.html"
    log:
        "workup/logs/multiqc.log"
    conda: 
        "envs/qc.yaml"
    shell: 
        "multiqc workup -o workup/qc"


# rule snpsplit_DNA:
#     input:
#         "workup/alignments/{sample}.DNA.all_bcs.bam"
#     output:
#         "workup/alignments/{sample}.DNA.all_bcs.allele_flagged.bam"
#     log:
#         "workup/logs/snpsplit.log"
#     conda:
#         "envs/snpsplit.yaml"
#     shell:
#         "SNPsplit --snp_file {snp_file} {input} &> {log}"




#from clusterflow pipeline
# we are currently using a very high penalty score for soft-clipping (--sp 1000,1000)
#because Hisat2 seems to soft-clip even when it should run in --end-to-end mode
# we are also filtering out unmapped reads (-F 4), or reads where the mate was unmapped (-F 8)
# we are also filtering non-primary alignments (-F 256)
#filter on mapq score of 20 (Skip alignments with MAPQ smaller than 20)
rule hisat2_align:
    input:
        fq="workup/fastqs/{sample}_R1.barcoded_rpm.fastq.gz"
    output:
        "workup/alignments/{sample}.hisat2.bam"
    threads: 8
    log:
        "workup/logs/{sample}.hisat2.log"
    shell:
        #--no-mixed \
        # --no-discordant $splices \
        "(hisat2 --sp 1000,1000 \
        -p 8 \
        -t \
        --phred33 \
        -x {hisat2_index} \
        -U {input.fq} | \
        samtools view -bSq 20 -F 4 -F 256 - > {output}) &> {log}"
#low alignment rate most likely due to not adaptor trimming and not allowing soft-clipping

rule bowtie2_align:
    input:
        fq="workup/fastqs/{sample}_R1.barcoded_dpm.fastq.gz"
    output:
        "workup/alignments/{sample}.bowtie2.bam"
    threads: 8
    log:
        "workup/logs/{sample}.bowtie2.log"
    shell:
        "(bowtie2 \
        -p 8 \
        -t \
        --phred33 \
        -x {bowtie2_index} \
        -U {input.fq} | \
        samtools view -bSq 20 -F 4 -F 256 - > {output}) &> {log}"


# unique mappers all have 255 MAPQ score
rule unique_mappers:
     input:
         "workup/alignments/{sample}.Aligned.sortedByCoord.out.bam"
     output:
         "workup/alignments/{sample}.Aligned.sortedByCoord.out.unique.bam"
     shell:
         "samtools view -b -q 255 -o {output} {input}"
