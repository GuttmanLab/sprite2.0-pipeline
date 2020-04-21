#!/usr/bin/env Rscript

if(suppressMessages(!require(argparse, quietly = T))){
  install.packages("argparse")
  suppressMessages(library(argparse, quietly = TRUE))
}

if(suppressMessages(!require(GenomicAlignments, quietly = T))){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("GenomicAlignments")
  suppressMessages(library(GenomicAlignments, quietly = TRUE))
}

if(suppressMessages(!require(rtracklayer, quietly = T))){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("rtracklayer")
  suppressMessages(library(rtracklayer, quietly = TRUE))
}

#parse bam file, overlap reads with gene annotation, quantify sense vs antisense

parser <- ArgumentParser(description='Estimate DNA contamination')
parser$add_argument('-b','--bam', metavar='bam', type="character", nargs='+', dest='bam',
                    help='BAM files for which to calculate % of reads on sense strand')
parser$add_argument('-g', '--gtf', dest='gtf', type='character',
                    help='GTF files path with transcript annotations')
parser$add_argument('-o', '--out', dest='out', type='character',
                    help='Output txt with % results')

args <- parser$parse_args()


#'Quantify the percentage of reads on sense strand
#'
#'Depending on library type would expect close to 100% or 0% for RNA-seq
#'
#'@param bam_path BAM file to be quantified
#'@param annotation_gr GTF genomicranges object
#'
#'@return percentage on sense strand
#'
sense_strand_percentage <- function(bam_path, annotation_gr){
  
  param <- ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE, isDuplicate=FALSE))
  aln <- GenomicAlignments::readGAlignments(bam_path, param = param)
  
  #To convert RNA-seq reads which may span exons is a little more difficult. 
  #Since a read can potentially span multiple exons,a single read may need to be converted to multiple ranges.
  #To solve this we can use the grglist() function to return a GRangesList with a separate GRanges for each read.
  grl <- grglist(aln,use.mcols = TRUE)

  #keep only main chromosomes
  gr_chrom <- keepStandardChromosomes(grl, pruning.mode="coarse")
  annotation_gr <- keepStandardChromosomes(annotation_gr, pruning.mode="coarse")
  annotation_gr <- diffloop::addchr(annotation_gr)

  
  #Query within subject (reads within genes)
  sense_ovlp <- findOverlaps(gr_chrom, annotation_gr, ignore.strand=FALSE, type="any")
  strand(annotation_gr) <- ifelse(strand(annotation_gr) == '+', '-', '+')
  anti_sense_ovlp <- findOverlaps(gr_chrom, annotation_gr, ignore.strand=FALSE, type="any")
  asl <- anti_sense_ovlp[!(queryHits(anti_sense_ovlp) %in% queryHits(sense_ovlp))]
  asl_o <- gr_chrom[queryHits(asl)]
  
  #Count each base
  # sense_len <- sum(width(unique(unlist(gr_chrom[queryHits(sense_ovlp[!duplicated(sense_ovlp)])]))))
  # anti_sense_len <- sum(width(unique(unlist(asl_o[!duplicated(asl_o)]))))
  
  #Remove hit duplicates (count each read only once) and PCR duplicates
  sense <- length(unique(unlist(gr_chrom[queryHits(sense_ovlp)])))
  #anything that doesn't already belong to a correct orientation gene
  anti_sense <- length(unique(unlist(asl_o)))
  
  perc_sense <- (sense/(sense+anti_sense))*100
  # ovlp_perc <- (sense_len/(sense_len+anti_sense_len))*100
  
  return(c(perc_sense, sense, anti_sense)) #ovlp_perc, sense_len, anti_sense_len
  
}

annotation_gr <- import(args$gtf)
gene_anno_gr <- annotation_gr[mcols(annotation_gr)$type == "gene"] #transcript
# chr16:11136592-11176393
# chr16 11144125-11144181
find_range <- GRanges(seqnames = "chr16", ranges = IRanges(start = 11144125, end = 11144181), strand = "+")
gene <- findOverlaps(find_range, gene_anno_gr, ignore.strand=FALSE)
reads <- findOverlaps(find_range, grl, ignore.strand=TRUE)
unlist(grl[subjectHits(reads)])




results <- list()
for(i in args$bam){
  file_basename <- basename(i)
  print(file_basename)
  results[[file_basename]] <- sense_strand_percentage(i, gene_anno_gr)
  
}

output <-plyr::ldply(results, rbind)
names(output) <- c("File", "Percentage_on_sense", "Sense", "Anti-sense")

write.table(output, args$out, sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#'Count number of overlaps within each gene
#'
#'@param bam_path Path to bam file
#'@param annotation_gr Intervals within which to count reads (needs to have gene_name)
#'
#'@return data.frame with annotation gene_name and counts
#'
count_within_gene_overlap <- function(bam_path, annotation_gr){
  
  param <- ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE))
  
  aln <- GenomicAlignments::readGAlignments("/mnt/data/2019_08_31_Merge_SPRITE_MegaSPRITE_HiSeq_Peter/workup/alignments/PB49.RNA.chr.bam", param = param)
  #To convert RNA-seq reads which may span exons is a little more difficult. 
  #Since a read can potentially span multiple exons,a single read may need to be converted to multiple ranges.
  # To solve this we can use the grglist() function to return a GRangesList with a separate GRanges for each read.
  grl <- grglist(aln,use.mcols = TRUE)
  # gr <- granges(aln, use.names=TRUE, use.mcols=TRUE)
  
  #keep only main chromosomes
  gr_chrom <- keepStandardChromosomes(grl, pruning.mode="coarse")
  annotation_gr <- keepStandardChromosomes(annotation_gr, pruning.mode="coarse")
  annotation_gr <- diffloop::addchr(annotation_gr)
  
  # sense_count <- countOverlaps(annotation_gr, gr_chrom, ignore.strand=FALSE, type="any")
  # strand(annotation_gr) <- ifelse(strand(annotation_gr) == '+', '-', '+')
  # anti_sense_count <- countOverlaps(annotation_gr, gr_chrom, ignore.strand=TRUE, type="any")
  # 
  #Query within subject (reads within genes)
  sense_ovlp <- findOverlaps(gr_chrom, annotation_gr)
  strand(annotation_gr) <- ifelse(strand(annotation_gr) == '+', '-', '+')
  anti_sense_ovlp <- findOverlaps(gr_chrom, annotation_gr)
  total <- findOverlaps(gr_chrom, annotation_gr, type = "within", ignore.strand=TRUE)
  asl <- anti_sense_ovlp[!(queryHits(anti_sense_ovlp) %in% queryHits(sense_ovlp))]
  asl_o <- gr_chrom[queryHits(asl)]
  # 
  # sense_hits <- duplicated(unlist(gr_chrom[queryHits(sense_ovlp)]))
  # anti_sense_hits <- duplicated(unlist(asl_o))
  # total_hits <- duplicated(unlist(gr_chrom[queryHits(total)]))
  sense_df <- data.frame(table(annotation_gr$gene_name[(subjectHits(sense_ovlp))]))
  anti_sense_df <- data.frame(table(annotation_gr$gene_name[(subjectHits(asl))]))
  total_df <- data.frame(table(annotation_gr$gene_name[(subjectHits(total))]))
  out_df <- merge(sense_df, anti_sense_df, by="Var1", all=TRUE)
  # annotation_gr$sense_count <- read_in_anno_count
  # annotation_gr$all_count <- read_in_anno_all_count

  
  return(out_df)
  
}


#'Quantify the percentage of reads on sense strand
#'
#'Depending on library type would expect close to 100% or 0% for RNA-seq
#'
#'@param bam_path BAM file to be quantified
#'@param annotation_gr GTF genomicranges object
#'
#'@return percentage on sense strand
#'
sense_strand_percentage_by_length <- function(bam_path, annotation_gr){
  param <- ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE))
  
  aln <- GenomicAlignments::readGAlignments(bam_path, param = param)
  
  grl <- grglist(aln,use.mcols = TRUE)
  
  #keep only main chromosomes
  gr_chrom <- keepStandardChromosomes(grl, pruning.mode="coarse")
  annotation_gr <- keepStandardChromosomes(annotation_gr, pruning.mode="coarse")
  annotation_gr <- diffloop::addchr(annotation_gr)
  
  
  #Query within subject (reads within genes)
  sense_ovlp <- findOverlaps(gr_chrom, annotation_gr, ignore.strand=FALSE, type="within")
  sense_isect <- pintersect(gr_chrom[queryHits(sense_ovlp)], annotation_gr[subjectHits(sense_ovlp)])
  sense_len <- sum(width(unique(unlist(sense_isect))))
  strand(annotation_gr) <- ifelse(strand(annotation_gr) == '+', '-', '+')
  anti_sense_ovlp <- findOverlaps(gr_chrom, annotation_gr, ignore.strand=FALSE, type="within")
  asl <- anti_sense_ovlp[!(queryHits(anti_sense_ovlp) %in% queryHits(sense_ovlp))]
  asl_o <- gr_chrom[queryHits(asl)]
  
  # sense_gr <- unique(unlist(gr_chrom[queryHits(sense_ovlp[!duplicated(sense_ovlp)])]))
  # anti_sense_gr <- unique(unlist(asl_o[!duplicated(asl_o)]))
  
  anti_sense_isect <- pintersect(gr_chrom[queryHits(asl)], annotation_gr[subjectHits(asl)])
  anti_sense_len <- sum(width(unique(unlist(anti_sense_isect))))
  
  perc_len <- (sense_len/(sense_len+anti_sense_len))*100
  return(c(perc_len, sense_len, anti_sense_len))
  
}
