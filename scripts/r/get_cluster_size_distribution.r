#!/usr/bin/env Rscript

#' @author Peter Chovanec, Noah Ollikainen, Charlotte A Lai

if(!require(ggplot2)){
    install.packages("ggplot2", repos='http://cran.us.r-project.org')
    library(ggplot2)
}

if(!require(ggpubr)){
    install.packages("ggpubr", repos='http://cran.us.r-project.org')
    library(ggpubr)
}

if(!require(reshape2)){
    install.packages("reshape2", repos='http://cran.us.r-project.org')
    library(reshape2)
}

#Data holder class
setClass(Class = "PlotData",
         slots=c(
             sizesLists="list",
             totalRnaCount="data.frame",
             totalDnaCount="data.frame",
             totalRRna="data.frame",
             countsLists="list",
             rpmOnlyClusters="data.frame",
             dpmOnlyClusters="data.frame",
             totalClusters="data.frame",
             rpmOnlyReads="data.frame",
             dpmOnlyReads="data.frame"),
         prototype = list(
             sizesLists=list(NULL),
             totalRnaCount=data.frame(NULL),
             totalDnaCount=data.frame(NULL),
             totalRRna=data.frame(NULL),
             countsLists=list(NULL),
             rpmOnlyClusters=data.frame(NULL),
             rpmOnlyClusters=data.frame(NULL),
             totalClusters=data.frame(NULL),
             rpmOnlyReads=data.frame(NULL),
             dpmOnlyReads=data.frame(NULL))
         )

setGeneric("sizesLists", function(object) standardGeneric("sizesLists"))
setMethod("sizesLists", "PlotData", function(object) object@sizesLists)

setGeneric("countsLists", function(object) standardGeneric("countsLists"))
setMethod("countsLists", "PlotData", function(object) object@countsLists)

setGeneric("totalRnaCount", function(object) standardGeneric("totalRnaCount"))
setMethod("totalRnaCount", "PlotData", function(object) object@totalRnaCount)

setGeneric("totalDnaCount", function(object) standardGeneric("totalDnaCount"))
setMethod("totalDnaCount", "PlotData", function(object) object@totalDnaCount)

setGeneric("totalRRna", function(object) standardGeneric("totalRRna"))
setMethod("totalRRna", "PlotData", function(object) object@totalRRna)


setGeneric("rpmOnlyClusters", function(object) standardGeneric("rpmOnlyClusters"))
setMethod("rpmOnlyClusters", "PlotData", function(object) object@rpmOnlyClusters)

setGeneric("dpmOnlyClusters", function(object) standardGeneric("dpmOnlyClusters"))
setMethod("dpmOnlyClusters", "PlotData", function(object) object@dpmOnlyClusters)

setGeneric("totalClusters", function(object) standardGeneric("totalClusters"))
setMethod("totalClusters", "PlotData", function(object) object@totalClusters)

setGeneric("rpmOnlyReads", function(object) standardGeneric("rpmOnlyReads"))
setMethod("rpmOnlyReads", "PlotData", function(object) object@rpmOnlyReads)

setGeneric("dpmOnlyReads", function(object) standardGeneric("dpmOnlyReads"))
setMethod("dpmOnlyReads", "PlotData", function(object) object@dpmOnlyReads)


SINGLETON        <- 1
FROM_2_TO_10     <- 2
FROM_11_TO_100   <- 3
FROM_101_TO_1000 <- 4
OVER_1000        <- 5

getCount <- function(x, n) x[[n]]

incrementCount <- function(x, cat, num) {
    x[[cat]] <-  x[[cat]] + num
    x
}

normalizeCount  <- function(x) {
    total <- sum(unlist(x))
    lapply(x, function(l) l / total)
}

catergoriseCount <- function(sizes, numReads){
    if      (numReads == 1)    sizes <- incrementCount(sizes, SINGLETON, numReads)
    else if (numReads <= 10)   sizes <- incrementCount(sizes, FROM_2_TO_10, numReads)
    else if (numReads <= 100)  sizes <- incrementCount(sizes, FROM_11_TO_100, numReads)
    else if (numReads <= 1000) sizes <- incrementCount(sizes, FROM_101_TO_1000, numReads)
    else                       sizes <- incrementCount(sizes, OVER_1000, numReads)
}

catergoriseRDpm <- function(sizes, numReads, numRDpm){
    if      (numReads == 1)    sizes <- incrementCount(sizes, SINGLETON, numRDpm)
    else if (numReads <= 10)   sizes <- incrementCount(sizes, FROM_2_TO_10, numRDpm)
    else if (numReads <= 100)  sizes <- incrementCount(sizes, FROM_11_TO_100, numRDpm)
    else if (numReads <= 1000) sizes <- incrementCount(sizes, FROM_101_TO_1000, numRDpm)
    else                       sizes <- incrementCount(sizes, OVER_1000, numRDpm)
}


makeDataFrame <- function(sizes, f, normalise=TRUE){
    if(normalise){
        sizes <- normalizeCount(sizes)   
    }
    df <- data.frame(unlist(sizes))
    colnames(df) <- c("count")
    fs <- c("SINGLETON", "FROM_2_TO_10", "FROM_11_TO_100", "FROM_101_TO_1000", "OVER_1000")
    df$category <- factor(fs, levels = c(fs[1], fs[2], fs[3], fs[4], fs[5]))
    df$filename <- basename(f)
    return(df)
}

#TODO add total DPM reads. total DPM only cluser reads, total RPM only cluster reads, total RPM and DPM cluster reads

processFile <- function(f) {
    sizes <- list(0, 0, 0, 0, 0)
    sizesDpm <- list(0, 0, 0, 0, 0)
    sizesRpm <- list(0, 0, 0, 0, 0)
    sizesDpmOverall <- list(0, 0, 0, 0, 0)
    sizesRpmOverall <- list(0, 0, 0, 0, 0)
    sizesDpmOnlyClusterReads <- 0
    sizesRpmOnlyClusterReads <- 0
    rpmOnlyClusters <- 0
    dpmOnlyClusters <- 0
    totalClusters <- 0
    rrnaCount <- 0

    tryCatch({
        if(endsWith(f, 'gz')){
            con <- gzcon(file(f,open="rb"))
        }else{
            con <- file(f, open = 'r')   
        }
        while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
            
            splitLine <- unlist(strsplit(oneLine, "\t"))
            numReads <- length(splitLine) - 1  # First column is barcode
            numDpm <- sum(grepl("DPM", splitLine))
            numRpm <- sum(grepl("RPM", splitLine))
            
            totalClusters <- totalClusters + 1
            if(numDpm == 0){
                sizesRpmOnlyClusterReads <- sizesRpmOnlyClusterReads + numRpm
                rpmOnlyClusters <- rpmOnlyClusters +1
            }else if(numRpm == 0){
                sizesDpmOnlyClusterReads <- sizesDpmOnlyClusterReads + numDpm
                dpmOnlyClusters <- dpmOnlyClusters +1
            }
            
            rrnaCount <- rrnaCount + sum(grepl(paste(c("rRNA", "5S", "28S", "45S", "5.8S", "18S", "4.5S", "ITS1", "ITS2"), 
                                                     collapse="|"), splitLine))
            
            sizes <- catergoriseCount(sizes, numReads)
            sizesDpm <- catergoriseCount(sizesDpm, numDpm)
            sizesRpm <- catergoriseCount(sizesRpm, numRpm)
            sizesDpmOverall <- catergoriseRDpm(sizesDpmOverall, numReads, numDpm)
            sizesRpmOverall <- catergoriseRDpm(sizesRpmOverall, numReads, numRpm)
        }
    }, finally = {
        close(con)
    })

    sizesNormList <- list(total=makeDataFrame(sizes, f), DpmCount=makeDataFrame(sizesDpm, f), 
                          RpmCount=makeDataFrame(sizesRpm, f), DpmOverall=makeDataFrame(sizesDpmOverall, f), 
                          RpmOverall=makeDataFrame(sizesRpmOverall, f))
    sizesList <- list(total=makeDataFrame(sizes, f, F), DpmCount=makeDataFrame(sizesDpm, f, F), 
                      RpmCount=makeDataFrame(sizesRpm, f, F), DpmOverall=makeDataFrame(sizesDpmOverall, f, F), 
                      RpmOverall=makeDataFrame(sizesRpmOverall, f, F))
    
    out <- new("PlotData",sizesLists=sizesNormList,
        totalRnaCount=data.frame(count=sum(sizesList[["RpmCount"]]["count"]), filename=basename(f)),
        totalDnaCount=data.frame(count=sum(sizesList[["DpmCount"]]["count"]), filename=basename(f)),
        totalRRna=data.frame(count=rrnaCount, filename=basename(f)),
        countsLists=sizesList,
        dpmOnlyClusters=data.frame(count=dpmOnlyClusters, filename=basename(f)),
        rpmOnlyClusters=data.frame(count=rpmOnlyClusters, filename=basename(f)),
        totalClusters=data.frame(count=totalClusters, filename=basename(f)),
        rpmOnlyReads=data.frame(count=sizesRpmOnlyClusterReads, filename=basename(f)),
        dpmOnlyReads=data.frame(count=sizesDpmOnlyClusterReads, filename=basename(f)))
    
    return(out)
}


combineDataFrames <- function(dfs){
    categories <- c("total", "DpmCount", "RpmCount", "DpmOverall", "RpmOverall")
    dfs_out <- list()
    for(cat in categories){
        dfs_out[[cat]] <- do.call("rbind",lapply(dfs,function(x) x[[cat]]))
    }
    return(dfs_out)
}


parseArgs <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 2) stop("get_cluster_sizes <directory> <pattern>")
    args
}

plotCounts <- function(df, dir_out, name_out) {
    p <- ggplot(df, aes(x = filename, y = count, fill = category)) +
         geom_bar(stat = "identity") +
         ylab("% of reads") +
         theme_pubr(legend="right") +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
               axis.title.x=element_blank()) +
         facet_grid(~origin)
    
    png(paste(dir_out, "/", name_out, "_cluster_sizes.png", sep=""))
    print(p)
    graphics.off()
    
    pdf(paste(dir_out, "/", name_out, "_cluster_sizes.pdf", sep=""))
    print(p)
    graphics.off()
}

plotSup <- function(df1, df2, data_names, dir_out, out_name){
    master_df <-merge(df1, df2, by="filename", all=TRUE)
    plot_df <- data.frame(master_df$filename, 
                         (master_df$count.x/(master_df$count.x+master_df$count.y))*100,
                         (master_df$count.y/(master_df$count.x+master_df$count.y))*100)
    names(plot_df) <- c("filename", data_names)
    plot_df <- reshape2::melt(plot_df)
    
    p <- ggplot(plot_df, aes(x = filename, y = value, fill = variable)) +
        geom_bar(stat = "identity") +
        ylab("% of reads") +
        theme_pubr(legend="right") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              axis.title.x=element_blank())
    
    png(paste(dir_out, paste("/", out_name, ".png", sep=""), sep=""))
    print(p)
    graphics.off()
    
    pdf(paste(dir_out, paste("/", out_name,".pdf", sep=""), sep=""))
    print(p)
    graphics.off()
    
}


main <- function() {
    args <- parseArgs()
    # files <- list.files(path = "/mnt/data/RNA_DNA_SPRITE/20200315_RNA_DNA_Mario",
    # pattern = ".clusters.gz", full.names = TRUE, recursive = FALSE)
    files <- list.files(path = args[1], pattern = args[2], full.names = TRUE, recursive = FALSE)
    
    #Make % stacked barplots
    dfs <- lapply(files, processFile)
    dfs_percentage <- combineDataFrames(lapply(dfs, function(x) sizesLists(x)))
    dfs_counts <- combineDataFrames(lapply(dfs, function(x) countsLists(x)))
    #add extra column for facet_grid
    label <- c("D/RPM % reads", "DPM % reads", "RPM % reads", 
               "DPM % of D/RPM", "RPM % of D/RPM")
    dfs_plot_list <- mapply(cbind, dfs_percentage, "origin"=label, SIMPLIFY=F)
    plot_df <- do.call("rbind", dfs_plot_list)
    plotCounts(plot_df, args[1], args[2])
    
    #Write out raw count tables
    dfs_count_list <- mapply(cbind, dfs_counts, "origin"=label, SIMPLIFY=F)
    count_df <- do.call("rbind", dfs_count_list)
    write.table(count_df, paste(args[1], "/", args[2], "_cluster_counts.tsv", sep=""), sep = '\t', 
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    #Make total RNA DNA rRNA plots
    rna_df <- do.call("rbind", lapply(dfs, function(x) totalRnaCount(x)))
    dna_df <- do.call("rbind", lapply(dfs, function(x) totalDnaCount(x)))
    plotSup(rna_df, dna_df, c("RNA", "DNA"), args[1], paste(args[2], "_rna_dna_proportions"))
    
    total_rrna <- do.call("rbind", lapply(dfs, function(x) totalRRna(x)))
    plotSup(total_rrna, rna_df, c("rRNA", "RNA"), args[1], paste(args[2],"_rrna_rna_proportions"))
    
    #Make table of other metrics
    # rna_df <- do.call("rbind", lapply(dfs, function(x) totalRnaCount(x)))
    # 
    # new("PlotData",sizesLists=sizesNormList,
    #     totalRnaCount=data.frame(count=sum(sizesList[["RpmCount"]]["count"]), filename=basename(f)),
    #     totalDnaCount=data.frame(count=sum(sizesList[["DpmCount"]]["count"]), filename=basename(f)),
    #     totalRRna=data.frame(count=rrnaCount, filename=basename(f)),
    #     countsLists=sizesList,
    #     dpmOnlyClusters=data.frame(count=dpmOnlyClusters, filename=basename(f)),
    #     rpmOnlyClusters=data.frame(count=rpmOnlyClusters, filename=basename(f)),
    #     totalClusters=data.frame(count=totalClusters, filename=basename(f)),
    #     rpmOnlyReads=data.frame(count=rpmOnlyReads, filename=basename(f)),
    #     dpmOnlyReads=data.frame(count=dpmOnlyReads, filename=basename(f)))
}

main()

