library(myvariant)
library(magrittr)
library(IRanges)
library(plyr)
library(httr)
library(jsonlite)
library(myvariant)
library(VariantAnnotation)
library(data.table)
source("~/.virtualenvs/op2/omics_pipe/omics_pipe/scripts/annotateIndels.R")

args <- commandArgs(TRUE)
vcfPath <- args[1]
somatic <- args[2]
filter <- args[3]


.collapse <- function (...) {
     paste(unlist(list(...)), sep = ",", collapse = ",")
 }

annotateVariantsFromVcf <- function(vcf.path, somatic=FALSE, do_filter=TRUE){
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
    ## MyVariant.info annotations
    vcf <- readVcf(vcf.path, genome="hg19")
    if (dim(vcf)[1] == 0){
	message("VCF is empty!")
	return()}
    snp <- vcf[isSNV(vcf)]
    hgvs <- formatHgvs(snp, "snp")
    annos <- getVariants(hgvs, fields=c("dbsnp.rsid", "cadd.consequence", 
                                        #"dbnsfp.aa.pos", "dbnsfp.aa.ref", "dbnsfp.aa.alt", 
                                        "snpeff.ann.hgvs_p",
                                       # "cadd.gene.prot.protpos", "cadd.oaa", "cadd.naa",
                                      "cadd.gene.genename",
                                       "cosmic.cosmic_id", "cosmic.tumor_site", "exac.af",
                                       "dbnsfp.1000gp3.af", "cadd.phred",
                                       "dbnsfp.polyphen2.hdiv.rankscore", "dbnsfp.polyphen2.hdiv.pred", 
                                       "dbnsfp.mutationtaster.converted_rankscore", "dbnsfp.mutationtaster.pred"
                                       ))
      dp <- geno(snp)$DP
      row.names(dp) <- NULL
      ad <- geno(snp)$AD
      row.names(ad) <- NULL
      coverage <- cbind(data.frame(ad), dp)
      if(somatic==FALSE){
        names(coverage) <- c("AD", "DP")
	}
      annos <- data.frame(cbind(annos, coverage))
    ## filter consequence
    if(do_filter){
      annos <- subset(annos, cadd.consequence %in% c("STOP_GAINED","STOP_LOST", 
                                                   "NON_SYNONYMOUS", "SPLICE_SITE", 
                                                   "CANONICAL_SPLICE", "REGULATORY"))
      annos <- arrange(data.frame(subset(annos, is.na(exac.af) | exac.af < 0.05)), -dbnsfp.mutationtaster.converted_rankscore)
    }
    setnames(annos, 
             old = c("query", "dbsnp.rsid", "cadd.consequence", "snpeff.ann",
                     "cosmic.cosmic_id", "cosmic.tumor_site", "exac.af", "dbnsfp.1000gp3.af", "cadd.phred",
                     "dbnsfp.polyphen2.hdiv.rankscore", "dbnsfp.polyphen2.hdiv.pred",
                     "dbnsfp.mutationtaster.converted_rankscore", "dbnsfp.mutationtaster.pred"), 
             new = c("Variant", "dbSNP rsid", "Consequence", "Amino Acid",
                     "COSMIC ID", "COSMIC Tumor Site", "ExAC AF", "G1000_AF", "CADD Score",
                     "Polyphen-2 Score", "Polyphen-2 Prediction", "MutationTaster Score", "MutationTaster Prediction"))
    names(annos)[names(annos) %in% c("cadd.gene.genename", "cadd.gene")] <- "Gene"
    annos <- DataFrame(annos) ##for some reason have to do this to eliminate the following columns
    annos[c("X_id", "notfound", "X_score", "cadd._license")] <- NULL
    annos <- lapply(annos, function(i) sapply(i, .collapse))
    ## merge annotations
    annotationsIntogen <- merge(annos, intogen, all.x=TRUE)
    annotationsPanCancer <- merge(annotationsIntogen, panCancer, all.x=TRUE)
    annotationsDruggability <- merge(annotationsPanCancer, druggability, all.x=TRUE, sort=TRUE)
    ## write file
    annotations <- data.frame(sapply(annotationsDruggability, as.character), stringsAsFactors = FALSE)
#    annotations <- rbind.fill(annotations, annotateIndels(vcf.path))
    annotationsFinal <- subset(annotations, !is.na(reported))
    annotationsFinal[is.na(annotationsFinal)] <- ""
    write.table(annotationsFinal, gsub(".vcf.gz", ".annotated.txt", vcf.path), sep="\t", row.names=FALSE, quote=FALSE)
    annotationsFinal
}

annotateVariantsFromVcf(vcfPath, somatic=somatic, do_filter=filter)
