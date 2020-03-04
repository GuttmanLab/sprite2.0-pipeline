
library(readr)

ligation_logs <- list.files("/mnt/data/RNA_DNA_SPRITE/2019_09_12_NovaSeq_MegaSPRITE_Peter/", "HL522*", full.names = TRUE)

#' Convert a list of ligation efficiency logs into a dataframe with full barcode %
#' 
#' @param list_of_files a list of file paths to ligation logs
#' @param num_barcodes how many barcodes constitute a full length barcodes
#' 
#' @return data.frame of full length barcode %
#' 
ligation_efficiency <- function(list_of_files, num_barcodes=5){
  
  results <- list()
  for(file in list_of_files){
    le <- readr::read_table2(file,col_names = FALSE)
    sample_name <- gsub(".ligation_efficiency.txt","", basename(file))
    perc_full <- gsub("\\(|\\)|%", "", le$X2[num_barcodes+1])
    results[[sample_name]] <- perc_full  
  }
  out <- plyr::ldply(results, rbind)  
  
  return(out)
}

ligation_df <- ligation_efficiency(ligation_logs, 5)




1. #/% of RNA DNA reads in cluster files. 
2. # of DNA reads pre and post masking. 
3. # of reads aligned for RNA and DNA. use mqc_bowtie2_se_plot_1.txt from multiqc file

  
alignment_logs <- function(){
  
  
}