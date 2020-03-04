#incomplete (can't read XT flags from bam file)

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


parser <- ArgumentParser(description='Estimate DNA contamination')
parser$add_argument('-b','--bam', metavar='bam', type="character",
                    help='BAM file for which to calculate DNA contamination score')
parser$add_argument('-g', '--gtf', dest='gtf', type='character',
                    help='sum the integers (default: find the max)')

args <- parser$parse_args()

# RNA_aln <- readGAlignments("D:/RNA_DNA_SPRITE/HL522DSXX.AACAATGG_TAGTTAGC.RNA.hisat2.mapq20.anno.bam")
# genes <- import("D:/genomes/GRCm38.p6/GRCm38.p6.annotation.gtf.gz")

RNA_aln <- readGAlignments("/mnt/data/RNA_DNA_SPRITE/HL522DSXX.AACAATGG_TAGTTAGC.RNA.hisat2.mapq20.anno.bam")
RNA_DNA_aln <- readGAlignments("/mnt/data/2019_08_31_Merge_SPRITE_MegaSPRITE_HiSeq_Peter/workup/alignments/PB49.RNA.hisat2.mapq20.anno.bam")
genes <- import("/mnt/data/genomes/GRCm38.p6/GRCm38.p6.annotation.gtf.gz")

# RNA_aln <- readGAlignments(args$bam) 
# genes <- import(args$gtf)

genes_anno <- genes[mcols(genes)$type == "gene"]
# 
RNA_gr <- granges(RNA_aln, use.names=TRUE, use.mcols=TRUE)
RNA_gr <- diffloop::rmchr(RNA_gr)
# #keep only main chr
chroms <- c(seq(1,19), "X", "Y")
RNA_gr_chrom <- RNA_gr[seqnames(RNA_gr) %in% chroms]
# 
# total_genome <- sum(seqlengths(seqinfo(RNA_gr_chrom)))
# 
# 
# intergenic_reads <- length(RNA_gr_chrom) - sum(countOverlaps(RNA_gr_chrom, genes_anno))
# 
# 
intergenic <- gaps(genes_anno)
intergenic_chrom <- intergenic[seqnames(intergenic) %in% chroms]
# total_length <- sum(width(intergenic_chrom))
# total_gene <- sum(width(genes_anno))
# 
# intergenic_reads_per_kb <- sum(countOverlaps(intergenic, RNA_gr_chrom))/(total_length/1000)
# gene_reads_per_kb <- sum(countOverlaps(genes_anno, RNA_gr_chrom))/(total_gene/1000)
# 
# percent_DNA <- (intergenic_reads_per_kb/gene_reads_per_kb)*100



#tile intergenic regions

intergenic_tile <- unlist(tile(intergenic_chrom, width=1000))
mcols(intergenic_tile)$count <- countOverlaps(intergenic_tile, RNA_gr_chrom, ignore.strand=TRUE)

#calculate density
mcols(intergenic_tile)$density <- mcols(intergenic_tile)$count/width(intergenic_tile)

dna_density_per_kb <- median(mcols(intergenic_tile)$density)

predicted_contamination <- dna_density_per_kb*(total_genome/1000)

cat("Predicted contamination of", args$bam, "is:", predicted_contamination, "\n")

test <- data.frame(intergenic_tile)
ggplot(test, aes(density)) + stat_ecdf(geom = "step")

# predictedContamination = (long)(dnaDensityPerKb[d]*(SeqMonkApplication.getInstance().dataCollection().genome().getTotalGenomeLength()/1000))


# table(strand(RNA_gr))
# 
# DNA_contamination <- RNA_aln[!RNA_aln %over% genes_anno,]
# RNA<- RNA_aln[RNA_aln %over% genes_anno,]
# 
# 
# length(DNA_contamination)/length(RNA_aln)







#'ggplot the coverage of a single chromosome for general qc
#'
#'@param gr genomicranges object with intervals to count (need seqlengths of chromosomes)
#'@param chrom the chromosome to plot
#'@param bin_size size of bins to count overlapping intervals over
#'
coverage_line_plot <- function(gr, chrom, bin_size=1000){
  
  kb_bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(gr),
                                       tilewidth=1000,
                                       cut.last.tile.in.chrom=TRUE)
  
  if("*" %in% strand(gr)){
    mcols(kb_bins)$count <- countOverlaps(kb_bins, gr)
    to_plot <- data.frame(kb_bins)
    
    p1 <- ggplot(to_plot[to_plot$seqnames == chrom,], aes(x=start, y=count)) +
      geom_line()
    
  } else {
    strand(kb_bins) <- "+"
    mcols(kb_bins)$count_pos <- countOverlaps(kb_bins, gr)
    strand(kb_bins) <- "-"
    mcols(kb_bins)$count_neg <- countOverlaps(kb_bins, gr)
    
    to_plot <- data.frame(kb_bins)
    
    pos <- to_plot[c("seqnames", "start", "count_pos")]
    neg <- to_plot[c("seqnames", "start", "count_neg")]
    names(pos)[3] <- "count"
    names(neg)[3] <- "count"
    pos$strand <- "+"
    neg$strand <- "-"
    
    plot <- rbind(pos, neg)
    
    p1 <- ggplot(plot[plot$seqnames == chrom,], aes(x=start, y=count, color=strand)) +
      facet_wrap(~strand) +
      geom_line()
    
  } 

  return(p1)
}

coverage_line_plot(RNA_gr_chrom, "X", 1000)



#' Cluster file parser for RNA DNA SPRITE
#' 
#' @param path File path to clusters file
#' 
#' @return Dataframe with first column being the barcode, second clusters
#' 
#' 
#' 
parse_clusters <- function(path, cores=8){
  
  clusters <- readr::read_table(path, col_names = FALSE)
  
  counts <- mclapply(clusters$X1, mc.cores=8, function(x){
    cl_row <- unlist(strsplit(x, "\t"))
    RNA <- cl_row[grepl("RPM", cl_row)]
    DNA <- cl_row[grepl("DPM", cl_row)]
    count <- length(cl_row[-1])
    DNA_count <- length(DNA)
    RNA_count <- length(RNA)
    return(c(cl_row[1], paste(DNA, collapse = "\t"), paste(RNA, collapse = "\t"), count, DNA_count, RNA_count))
  })
  
  out <- plyr::ldply(counts, rbind)
  names(out) <- c("Barcode","DNA", "RNA", "cluster_size", "DNA_count", "RNA_count")
  return(out)
}


pb49_clusters <- parse_clusters("/mnt/data/2019_08_31_Merge_SPRITE_MegaSPRITE_HiSeq_Peter/workup/clusters/PB49.clusters")


library(stringr)
library(data.table)


gene_strand_lookup <- data.table(names = genes$gene_name, strand = as.character(strand(genes)), key = "names")
gene_strand_lookup <- unique(gene_strand_lookup)


#'Get strand of gene within cluster read
#'
#'@param cluster_read a single read form a cluster (should only have a single gene)
#'@param strand_lookup data.table with gene_name and strand, with gene_name being the key
#'
#'@return + or -
#'
find_gene_strand <- function(cluster_read, strand_lookup){
  exon_intron <- stringr::str_extract(cluster_read, "\\w+(\\.intron|\\.exon)")
  cl_genes <- unique(unlist(strsplit(exon_intron, "\\.intron|\\.exon")))
  stopifnot(length(cl_genes) == 1)
  gene_strand <- gene_strand_lookup[cl_genes]$strand
  return(gene_strand)
  
}

find_gene_strand("RPM[Adam32.intron;Unassigned_NoFeatures.none;MTC,LTR,ERVL-MaLR,MTC_dup21836.repeat]_chr8:24927495-24927580", gene_strand_lookup)





