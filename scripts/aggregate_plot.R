
library(EnrichedHeatmap)
library(readr)
library(ggplot2)
library(rtracklayer)

rna_7sk <- rtracklayer::import("/mnt/data/RNA_DNA_SPRITE/7SK.100bps.bedgraph", format = "bedGraph")

ES_log2_FPKM <- read_delim("/mnt/data/RNA_DNA_SPRITE/ES_log2_FPKM.txt", 
                           "\t", escape_double = FALSE, col_types = cols(Chr = col_character()), 
                           trim_ws = TRUE)



#Highly expressed (red, fragments per kilobase of transcript per million [FPKM] > 10), moderately expressed (gray FPKM = 2–10), or inactive (blue, FPKM = 0–2) genes are indicated.

qplot(ES_log2_FPKM$`Mean doxminusAligned.sortedByCoord.out.bam`)
qplot(ES_log2_FPKM$`Mean total.gsnap.seAligned.sortedByCoord.out.bam`)



highly_active_genes <- ES_log2_FPKM[ES_log2_FPKM$`Mean doxminusAligned.sortedByCoord.out.bam` > 5,]
moderatly_active_genes <- ES_log2_FPKM[ES_log2_FPKM$`Mean doxminusAligned.sortedByCoord.out.bam` < 5 & ES_log2_FPKM$`Mean doxminusAligned.sortedByCoord.out.bam` > 0,]
inactive_genes <- ES_log2_FPKM[ES_log2_FPKM$`Mean doxminusAligned.sortedByCoord.out.bam` < 0,]

highly_active_genes$type <- "High"
moderatly_active_genes$type <- "Moderate"
inactive_genes$type <- "Inactive"


hi_active_gr <- makeGRangesFromDataFrame(highly_active_genes[c("Chr", "Start", "End", "Strand", "Feature", "type")], keep.extra.columns=TRUE)
hi_active_gr_chr <- diffloop::addchr(hi_active_gr)
mo_active_gr <- makeGRangesFromDataFrame(moderatly_active_genes[c("Chr", "Start", "End", "Strand", "Feature", "type")], keep.extra.columns=TRUE)
mo_active_gr_chr <- diffloop::addchr(mo_active_gr)
in_active_gr <- makeGRangesFromDataFrame(inactive_genes[c("Chr", "Start", "End", "Strand", "Feature", "type")], keep.extra.columns=TRUE)
in_active_gr_chr <- diffloop::addchr(in_active_gr)

export.bed(hi_active_gr_chr, "/mnt/data/RNA_DNA_SPRITE/hi_active_gr_chr.bed")
export.bed(mo_active_gr_chr, "/mnt/data/RNA_DNA_SPRITE/mo_active_gr_chr.bed")
export.bed(in_active_gr_chr, "/mnt/data/RNA_DNA_SPRITE/in_active_gr_chr.bed")

#mm9
ES_log2_FPKM_mm9 <- read_delim("/mnt/data/RNA_DNA_SPRITE/ES_fpkm_mm9.txt", 
                           "\t", escape_double = FALSE, col_types = cols(Chr = col_character()), 
                           trim_ws = TRUE)

qplot(ES_log2_FPKM_mm9$`Mean total.gsnap.sorted.bam`)



highly_active_genes_mm9 <- ES_log2_FPKM_mm9[ES_log2_FPKM_mm9$`Mean total.gsnap.sorted.bam` > 5,]
moderatly_active_genes_mm9 <- ES_log2_FPKM_mm9[ES_log2_FPKM_mm9$`Mean total.gsnap.sorted.bam` < 5 & ES_log2_FPKM_mm9$`Mean total.gsnap.sorted.bam` > 0,]
inactive_genes_mm9 <- ES_log2_FPKM_mm9[ES_log2_FPKM_mm9$`Mean total.gsnap.sorted.bam` < 0,]

highly_active_genes_mm9$type <- "High"
moderatly_active_genes_mm9$type <- "Moderate"
inactive_genes_mm9$type <- "Inactive"


hi_active_gr <- makeGRangesFromDataFrame(highly_active_genes_mm9[c("Chr", "Start", "End", "Strand", "Feature", "type")], keep.extra.columns=TRUE)
hi_active_gr_chr <- diffloop::addchr(hi_active_gr)
mo_active_gr <- makeGRangesFromDataFrame(moderatly_active_genes_mm9[c("Chr", "Start", "End", "Strand", "Feature", "type")], keep.extra.columns=TRUE)
mo_active_gr_chr <- diffloop::addchr(mo_active_gr)
in_active_gr <- makeGRangesFromDataFrame(inactive_genes_mm9[c("Chr", "Start", "End", "Strand", "Feature", "type")], keep.extra.columns=TRUE)
in_active_gr_chr <- diffloop::addchr(in_active_gr)

export.bed(hi_active_gr_chr, "/mnt/data/RNA_DNA_SPRITE/hi_active_gr_chr_mm9.bed")
export.bed(mo_active_gr_chr, "/mnt/data/RNA_DNA_SPRITE/mo_active_gr_chr_mm9.bed")
export.bed(in_active_gr_chr, "/mnt/data/RNA_DNA_SPRITE/in_active_gr_chr_mm9.bed")





all_gr <- c(hi_active_gr_chr, mo_active_gr_chr, in_active_gr_chr)

partition <- mcols(all_gr)$type
# naive_is <- naive_is[!(is.na(naive_is$score)|is.infinite(naive_is$score))]

#Active genes
# primed_is <- rtracklayer::import("/media/chovanec/My_Passport/CHiC_naive_primed/washu_tracks/primed_hESC-1_insulation_sorted.bedGraph.gz", format = "bedGraph")
# primed_is <- primed_is[!(is.na(primed_is$score)|is.infinite(primed_is$score))]
# 
# 
# naive_tads
# 

mat_all <- normalizeToMatrix(rna_7sk, all_gr, value_column = "score", mean_mode = "absolute",
                               extend = 1000, w = 100, background = 0, smooth = TRUE)

nm_col_fun <- circlize::colorRamp2(quantile(mat_all, c(0, 0.99)), c("white", "red"))
# nm_col_fun_mo <- circlize::colorRamp2(quantile(mat_ma, c(0, 0.99)), c("white", "red"))
# nm_col_fun_ia <- circlize::colorRamp2(quantile(mat_ia, c(0, 0.99)), c("white", "red"))

lgd <- Legend(at = c("High", "Moderate", "Inactive"), title = "Expresion", 
             type = "lines", legend_gp = gpar(col = 2:4))

ht_list <- EnrichedHeatmap(mat_all, col = nm_col_fun, name = "7sk", axis_name_rot = 90,
                           column_title = "7sk", split=partition,
                           top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4))))
draw(ht_list, main_heatmap = "7sk", annotation_legend_list = list(lgd))

