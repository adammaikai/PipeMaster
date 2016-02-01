# original author: Tobias Meissner
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(edgeR)
library(dplyr)
library(SPIA)
library(graphite)
library(org.Hs.eg.db)
library(ggplot2)
library(corrplot)
library(EDASeq)
library(Nof1RNASeq)
################################################################################
# get command line arguments
################################################################################
args <- commandArgs(trailingOnly = TRUE)

patientData <- args[1] #'~/AWS/storage/pats/rna/CCD018/CCD018.ReadsPerGene.out.tab' 
tumorType <- args[2] #'BRCA'
libType <- args[3] #'stranded-reverse' # options: unstranded, stranded, stranded-reverse
refPath <- args[4] #'~/AWS/s3/averaprojects/CRC'
outDir <- args[5] #'~/Downloads/temp/'

if(dir.exists(outDir)) {
  print('Output directory exists')
  } else {
    dir.create(outDir)
}

################################################################################
# load the reference cohort
################################################################################
source('src/loadRef.R')
allTumorTypes <- c('BRCA', 'OV', 'LAML',
                   'FALLOPIAN')
if (tumorType %in% allTumorTypes) {
  loadRef(tumorType, refPath)
} else {
  print('Tumor Type is not supported')
  stop()
}

vsdRefMat <- assay(reference$reference_vsd)

################################################################################
# read in patient data
################################################################################
source('src/patientDataIn.R')
df <- patData(patientData)

################################################################################
# diff expression against healthy controls and/or tcga matched normal
################################################################################
source('src/diff.R')
diffExpr <- diff(tumorType, vsdRefMat, df)

res1x <- as.data.frame(diffExpr$res1)
res1x <- cbind(id=rownames(res1x), res1x) #,@call=patientREC$Pcalls)
res1Sig <- arrange(res1x[which(res1x$padj < 0.05), ], padj)
write.table(res1Sig, paste0(outDir, 'patient_vs_normal.csv'), sep='\t', row.names = FALSE)

if(any(reference$group=='MNORMAL')) {
  res2x <- as.data.frame(diffExpr$res2)
  res2x <- cbind(id=rownames(res2x), res2x) #,@call=patientREC$Pcalls)
  res2Sig <- arrange(res2x[which(res2x$padj < 0.05), ], padj)
  write.table(res2Sig, paste0(outDir, 'patient_vs_mnormal.csv'), sep='\t', row.names = FALSE)  
}

################################################################################
# viz
################################################################################
source('src/viz.R')
palette(rainbow(4))

if(any(reference$group=='MNORMAL')) {
  pdf(paste0(outDir, 'qc.pdf'))
  # add sample to reference
  sRef <- samleToRef(df, 
                     reference, 
                     loggeomeansRef, 
                     vsdRefMat, 
                     diffExpr$selNormal, 
                     diffExpr$selMNormal, 
                     diffExpr$selTumor
  )
} else {
  pdf(paste0(outDir, 'qc.pdf'))
  # add sample to reference
  sRef <- samleToRef(df, 
                     reference, 
                     loggeomeansRef, 
                     vsdRefMat, 
                     diffExpr$selNormal, 
                     NULL,
                     diffExpr$selTumor
  )
}  

## library size
my.libplot(diffExpr$dds, col=factor(sRef$des), legend=sRef$des)
  
## transformed data, boxplot
my.vsdplot(sRef$vsdMat, col=factor(sRef$des), sRef$des)
  
## relative log expression
plotRLE(sRef$vsdMat, 
        col=factor(sRef$des), 
        outline=FALSE, 
        las=3, 
        ylim=c(-.2, .2), 
        ylab="Relative Log Expression", 
        cex.axis=1, 
        cex.lab=1,
        names=sRef$des
)
  
## pca
my.pca(sRef$vsdMat, sRef$des)
  
## sample-to-sample dist
my.ssDist(sRef$vsdMat, sRef$des)
dev.off()




