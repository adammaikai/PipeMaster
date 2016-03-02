library(DNAcopy)

args <- commandArgs(TRUE)
prefix <- args[1] #/data/storage/patients/cnv
sample <- args[2] #sampleID

cn <- read.table(paste0(prefix, "/", sample, "/", sample, "-cytoSNP850k"), header=F, sep=" ", stringsAsFactors = F)

CNA.object <-CNA( genomdat = cn[,4], chrom = cn[,2], maploc = cn[,3], data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, verbose=0, min.width=2)
seg.pvalue <- segments.p(segs, ngrid=100, tol=1e-6, alpha=0.05, search.range=100, nperm=1000)

write.table (seg.pvalue, paste0(prefix, "/", sample, "/", sample, ".cytoSNP850k.seg"), row.names=F, col.names=F, quote=F, sep="\t")
