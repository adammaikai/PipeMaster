library(DNAcopy)

args <- commandArgs(TRUE)
prefix <- args[1] ##/data/storage/patients/vcf
sample <- args[2] 

###cn <- read.table("/data/storage/patients/vcf/CCD025/CCD025.varscan_cnv.copynumber",header=T, stringsAsFactors = F)
cn <- read.table(paste0(prefix, sample, "/", sample, ".varscan_cnv.copynumber"), header=T, sep="\t", stringsAsFactors = F)

CNA.object <-CNA( genomdat = cn[,7], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, verbose=0, min.width=2)
seg.pvalue <- segments.p(segs, ngrid=100, tol=1e-6, alpha=0.05, search.range=100, nperm=1000)

###write.table (seg.pvalue, file="/data/storage/patients/vcf/CCD025/CCD025_varscan_cnv.copynumber.outfile", row.names=F, col.names=F, quote=F, sep="\t")
write.table (seg.pvalue, paste0(prefix, sample, "/", sample, ".varscan_cnv.copynumber.seg"), row.names=F, col.names=F, quote=F, sep="\t")
