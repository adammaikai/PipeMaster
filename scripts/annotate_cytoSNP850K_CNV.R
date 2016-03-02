library(plyr)
library(mygene)

args <- commandArgs(TRUE)
prefix = args[1] #/data/storage/patients/cnv
sample = args[2] #sampleID

.collapse <- function (...) 
{
    paste(unlist(list(...)), sep = ",", collapse = ",")
}

act <- read.csv("/data/database/druggability/cancer_genes.txt", sep="\t", stringsAsFactors = F)
## annotate Varscan CNV output, intersect with Actionable Genes curated by Tobias
annotate_cytoSNP850k_CNV <- function(sample){
  cnv <- read.csv(paste0(prefix, "/", sample, "/", sample, ".cytoSNP850k.events.tsv"), sep = "\t", stringsAsFactors = FALSE)
  cnv$query <- paste0("hg19.",cnv$chrom, ":", cnv$chr_start, "-", cnv$chr_stop)
  cnv$symbol <- sapply(cnv$query, function(i) .collapse(query(i, species="human")$hits$symbol[[1]]))
  final <- arrange(subset(cnv, symbol %in% act$Symbol), -abs(seg_mean))
  write.table(final, file=paste0(prefix, "/", sample, "/", sample, "_cytoSNP850k_cnv_annotations.txt"), sep="\t", row.names=FALSE, quote=FALSE)
  }
annotate_cytoSNP850k_CNV(sample)
