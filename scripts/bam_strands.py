import gtf


gtf_anno = gtf.parse_gtf("/mnt/data/genomes/GRCm38.p6/GRCm38.p6.annotation.gtf.gz", "gene")


#parse bam with RNA or DNA, find overlapping reads with genes based on strand (using gene names)