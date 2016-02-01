library(VariantAnnotation)
library(myvariant)
library(mygene)
library(data.table)

.collapse <- function (...) {
     paste(unlist(list(...)), sep = ",", collapse = ",")
 }

annotateIndels <- function(vcf.path, somatic=FALSE){
  ## intogen file
  intogen <- read.csv("/data/database/druggability/Mutational_drivers_per_tumor_type.tsv", sep="\t", comment.char="#")
  names(intogen) <- c("Gene", "tumor_type")
  intogen <- aggregate(tumor_type ~ Gene, intogen, .collapse)

  ## panCancer file
  panCancer <- read.csv("/data/database/druggability/NanostringPanCancerGeneList.csv", sep=",", comment.char="#")
  names(panCancer) <- c("Gene", "reported")
  panCancer <- subset(aggregate(reported ~ Gene, panCancer, .collapse), !grepl("reference", reported) & !grepl("control", reported))

  ## intogen druggability
  druggability <- read.csv("/data/database/druggability/Protein_Drug_Interactions.tsv", sep="\t", comment.char="#")
  druggability <- druggability[, c("geneHGNCsymbol", "Drug_name")]
  druggability <- aggregate(Drug_name ~ geneHGNCsymbol, druggability, .collapse)
  names(druggability) <- c("Gene", "drug_name")
  
  vcf.object <- readVcf(vcf.path, genome="hg19")
  indel <- vcf.object[isIndel(vcf.object)]
  if (dim(indel)[1] == 0){
    return(data.frame())
  }
  loc <- paste("hg19.", seqnames(indel), ":", start(indel), "-", end(indel), sep="")
  Variant <- formatHgvs(indel, variant_type=c("insertion", "deletion"))
  hits <- lapply(loc, function(i) query(i)$hits$symbol)
  
  dp <- data.frame(geno(indel)$DP)
  row.names(dp) <- NULL
  # Allelic Depth
  ad <- data.frame(geno(indel)$AD)
  row.names(ad) <- NULL
  if (grepl("varscan", vcf.path)){
    setnames(dp, old=c(names(dp)[grepl("TUMOR", names(dp))], names(dp)[grepl("NORMAL", names(dp))]), 
             new=c("TUMOR.DP", "NORMAL.DP"))
    setnames(ad, old=c(names(ad)[grepl("TUMOR", names(ad))], names(ad)[grepl("NORMAL", names(ad))]), 
             new=c("TUMOR.AD", "NORMAL.AD"))
    coverage <- cbind(ad, dp)
    coverage$NORMAL.ALT.AF <- round(coverage$NORMAL.AD/coverage$NORMAL.DP, 2)
    coverage$TUMOR.ALT.AF <- round(coverage$TUMOR.AD/coverage$TUMOR.DP, 2)
    coverage$VariantCaller <- "Varscan2 Somatic"
  }
  if (grepl("mutect", vcf.path)){
    setnames(dp, old=c(names(dp)[grepl(".T", names(dp))], names(dp)[grepl(".B", names(dp))]), 
             new=c("TUMOR.DP", "NORMAL.DP"))
    setnames(ad, old=c(names(ad)[grepl(".T", names(ad))], names(ad)[grepl(".B", names(ad))]), 
             new=c("TUMOR.AD", "NORMAL.AD"))
    coverage <- cbind(ad, dp)
    naf <- lapply(coverage$NORMAL.AD, function(i) as.vector(i[[2]]/(i[[1]] + i[[2]])))
    names(naf) <- NULL
    taf <- lapply(coverage$TUMOR.AD, function(i) as.vector(i[[2]]/(i[[1]] + i[[2]])))
    names(taf) <- NULL
    coverage$NORMAL.ALT.AF <- lapply(naf, round, 2)
    coverage$TUMOR.ALT.AF <- lapply(taf, round, 2)
    coverage$VariantCaller <- "MuTect"
  }
  annotations <- DataFrame(Variant=Variant, Gene=sapply(hits, .collapse))
  annotations <- cbind(annotations, coverage)
  annotationsIntogen <- merge(annotations, intogen, all.x=TRUE)
  annotationsPanCancer <- merge(annotationsIntogen, panCancer, all.x=TRUE)
  annotationsDruggability <- merge(annotationsPanCancer, druggability, all.x=TRUE, sort=TRUE)
  annotations <- data.frame(sapply(annotationsDruggability, as.character), stringsAsFactors = FALSE)
  annotations
  data.frame(subset(annotations, !is.na(tumor_type)))
}
