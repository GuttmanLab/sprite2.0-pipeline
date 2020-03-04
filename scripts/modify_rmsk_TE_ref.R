
#NOT USED ----------------------------------------------------------------------

library(rtracklayer)

rmsk_te <- import.gff("/mnt/data/genomes/GRCm38.p6/repeats/GRCm38_rmsk_TE.gtf.gz")

mcols(rmsk_te)$all_id <- apply(mcols(rmsk_te)[c("gene_id", "transcript_id", "family_id", "class_id")], 1, function(x) paste(x, collapse=","))

export.gff(rmsk_te, "/mnt/data/genomes/GRCm38.p6/repeats/GRCm38_rmsk_TE_all_id.gtf.gz")


rmsk_te_h <- import.gff("/mnt/data/genomes/GRCh38/GRCh38_rmsk_TE.gtf.gz")

mcols(rmsk_te_h)$all_id <- apply(mcols(rmsk_te_h)[c("gene_id", "transcript_id", "family_id", "class_id")], 1, function(x) paste(x, collapse=","))

export.gff(rmsk_te, "/mnt/data/genomes/GRCh38/GRCh38_rmsk_TE_all_id.gtf.gz")


#-------------------------------------------------------------------------------

mm10_anno <- import.gff("/mnt/data/genomes/GRCm38.p6/GRCm38.p6.annotation.gtf")
  
data.frame(mm10_anno[mcols(mm10_anno)$gene_name == "Nanog"])
  


intronsByTranscript
library(GenomicFeatures)
library(rtracklayer)

gtf <- makeTxDbFromGFF("/mnt/data/genomes/GRCm38.p6/GRCm38.p6.annotation.gtf.gz", format = "gtf")
introns <- intronsByTranscript(gtf)
export.gff(introns, "intergenic.gff", format="gff")

#introns by gene https://www.biostars.org/p/165226/
exons <- exonsBy(gtf, by="gene")

#make introns
exons <- reduce(exons)
exons <- exons[sapply(exons, length) > 1]

introns <- lapply(exons, function(x) {
  #Make a "gene" GRange object
  gr = GRanges(seqnames=seqnames(x)[1], ranges=IRanges(start=min(start(x)),
                                                       end=max(end(x))), 
               strand=strand(x)[1])
  db = disjoin(c(x, gr))
  ints = db[countOverlaps(db, x) == 0]
  #Add an ID
  if(as.character(strand(ints)[1]) == "-") {
    ints$exon_id = c(length(ints):1)
  } else {
    ints$exon_id = c(1:length(ints))
  }
  ints
})
introns <- GRangesList(introns)



gtf = makeTxDbFromGFF("mygtf.gtf", format = "gtf")

transcript = exonsBy(gtf, "tx")
intergenic = gaps(unlist(range(transcript)))
export.gff(intergenic, "intergenic.gff", format="gff")


exon_intron <- import.gff("/mnt/data/genomes/GRCm38.p6/GRCm38.p6.exon_intron.gtf")
