library(readr)
library(rtracklayer)
library(GenomicRanges)

#' Convert USCS rmsk txt file to genomicRanges object
#' 
#' UCSC uses 0-based coordinate system
#' 
#' @param rmsk_txt path to UCSC rmsk txt file
#' @param milliDev keep repeats younger or equal to (140)
#' 
#' @return GenomicRanges object
#' 
rmsk2gr <- function(rmsk_txt, milliDev=140){
  hg38_rmsk_txt <- read_delim(rmsk_txt, 
                              "\t", escape_double = FALSE, col_names = FALSE, 
                              trim_ws = TRUE)
  
  col_names <- c("bin",	"swScore",	"milliDiv",	"milliDel",	"milliIns",	"genoName",	
                 "genoStart",	"genoEnd",	"genoLeft",	"strand",	"repName",	"repClass",	
                 "repFamily",	"repStart",	"repEnd",	"repLeft",	"id")
  
  names(hg38_rmsk_txt) <- col_names
  
  rmsk_140milli <- hg38_rmsk_txt[hg38_rmsk_txt$milliDiv <= 140,]
  
  rmsk_140milli$name <- paste(rmsk_140milli$repName, rmsk_140milli$repClass, rmsk_140milli$repFamily, sep="|")
  bed_out <- rmsk_140milli[c("genoName",	"genoStart",	"genoEnd", "milliDiv", "name", "strand")]
  
  rmsk_gr <- makeGRangesFromDataFrame(bed_out, keep.extra.columns = TRUE, seqnames.field = ("genoName"), 
                                      starts.in.df.are.0based=TRUE)
  
  return(rmsk_gr)
}


hg38_rmsk_gr <- rmsk2gr("/mnt/data/genomes/GRCh38/repeats/hg38_rmsk.txt.gz")
mm10_rmsk_gr <- rmsk2gr("/mnt/data/genomes/GRCm38.p6/repeats/mm10_rmsk.txt.gz")
# write.table(bed_out, "/mnt/data/genomes/GRCh38/hg38_rmsk_140milliDev.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# Make sure blacklisted regions do not overlap ---------------------------------

#Previous blacklist had overlapping coordinates, just in case, perform merging
blacklist_mm10 <- rtracklayer::import.bed("/mnt/data/genomes/GRCm38.p6/mm10-blacklist.v2.bed.gz")
blacklist_hg38 <- rtracklayer::import.bed("/mnt/data/genomes/GRCh38/hg38-blacklist.v2.bed.gz")

blacklist_mm10_r <- reduce(blacklist_mm10)
mcols(blacklist_mm10_r)$name <- mcols(blacklist_mm10)[match(blacklist_mm10, blacklist_mm10_r), "name"]


blacklist_hg38_r <- reduce(blacklist_hg38)
mcols(blacklist_hg38_r)$name <- mcols(blacklist_hg38)[match(blacklist_hg38, blacklist_hg38_r), "name"]


blacklist_hg38_rmsk_milliDev140 <- c(hg38_rmsk_gr, blacklist_hg38_r)
blacklist_mm10_rmsk_milliDev140 <- c(mm10_rmsk_gr, blacklist_mm10_r)

rtracklayer::export.bed(blacklist_hg38_rmsk_milliDev140, "/mnt/data/genomes/GRCh38/hg38_blacklist_rmsk.milliDivLessThan140.bed")
rtracklayer::export.bed(blacklist_mm10_rmsk_milliDev140, "/mnt/data/genomes/GRCm38.p6/mm10_blacklist_rmsk.milliDivLessThan140.bed")

