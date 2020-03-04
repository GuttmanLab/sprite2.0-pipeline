library(tximport)
library(rhdf5)
library(rtracklayer)
library(diffloop)
library(dplyr)

gtf <- import.gff("/mnt/data/genomes/GRCm38.p6/GRCm38.p6.annotation.gtf")
tx2gene <- mcols(gtf)[c("transcript_id", "gene_name")]

gtf_gene <- gtf[gtf$type == "gene"]


files <- "/mnt/data/quant_DOXminus/abundance.h5"
txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion=TRUE)

tpm <- data.frame(txi.kallisto$abundance)
tpm$gene_name <- row.names(tpm)

tpm_nz <- tpm[tpm$txi.kallisto.abundance > 0,]

hist(log2(tpm_nz$txi.kallisto.abundance), breaks = 50)
abline(v=5,col="red")
abline(v=2,col="red")

tpm_high <- tpm_nz[log2(tpm_nz$txi.kallisto.abundance) >= 5,]
tpm_medium <- tpm_nz[log2(tpm_nz$txi.kallisto.abundance) >= 2 & log2(tpm_nz$txi.kallisto.abundance) < 5,]
tpm_low <- tpm_nz[log2(tpm_nz$txi.kallisto.abundance) <= 2,]

#'Convert gene tpm  list to genomeRanges objects with gene position
#'
#'@param tpm data.frame with gene_name and tpm columns
#'@param gtf_gr GenomicRanges object of annotation GTF
#'
gene2bed <- function(tpm, gtf_gr){
  
  
  gtf_gr_sub <- gtf_gr[gtf_gr$gene_name %in% tpm$gene_name]
  
  tpm_ord <- tpm %>% slice(match(gtf_gr_sub$gene_name, gene_name))
  
  #add tpm values to genomicRanges object
  gtf_gr_sub$tpm <- tpm_ord$txi.kallisto.abundance
  gr_chr <- diffloop::addchr(gtf_gr_sub)
  mcols(gr_chr) <- mcols(gr_chr)[c("gene_name","gene_biotype", "tpm")]
  
  return(gr_chr)
}

export.bed(gene2bed(tpm_high, gtf_gene), "/mnt/data/RNA_DNA_SPRITE/hi_active_kallisto_mm10.bed")
export.bed(gene2bed(tpm_medium, gtf_gene), "/mnt/data/RNA_DNA_SPRITE/med_active_kallisto_mm10.bed")
export.bed(gene2bed(tpm_low, gtf_gene), "/mnt/data/RNA_DNA_SPRITE/low_active_kallisto_mm10.bed")

write.table(data.frame(gene2bed(tpm_nz, gtf_gene)), "/mnt/data/RNA_DNA_SPRITE/DOXminus_kallisto_mm10.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
