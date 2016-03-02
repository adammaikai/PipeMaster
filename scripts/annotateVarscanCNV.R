library(plyr)
library(mygene)

args <- commandArgs(TRUE)
prefix = args[1] ###/data/storage/patients
sample = args[2]

.collapse <- function (...) 
{
    paste(unlist(list(...)), sep = ",", collapse = ",")
}

act <- read.csv("/data/database/druggability/actionable_avera.tsv", sep="\t", stringsAsFactors = F)
## annotate Varscan CNV output, intersect with Actionable Genes curated by Tobias
annotateVarscanCNV <- function(sample){
  cnv <- read.csv(paste0(prefix, "/vcf/", sample, "/", sample, ".varscan_cnv.events.tsv"), sep = "\t", stringsAsFactors = FALSE)
  cnv$query <- paste0(cnv$chrom, ":", cnv$chr_start, "-", cnv$chr_stop)
  cnv$symbol <- sapply(cnv$query, function(i) .collapse(query(i, species="human")$hits$symbol[[1]]))
  final <- arrange(subset(cnv, symbol %in% act$Symbol), -abs(seg_mean))
  write.table(final, file=paste0(prefix, "/cnv/", sample, "/", sample, "_varscan_cnv_annotations.txt"), sep="\t", row.names=FALSE, quote=FALSE)
  }
annotateVarscanCNV(sample)
