## Script to process data from the Nanostring CNV platform and call copy numbers for tumor-normal samples.

library(plyr)

args <- commandArgs(TRUE)
prefix = args[1] #/data/storage/patients/cnv
sample = args[2] #sampleID

tfile <- paste0(prefix, "/", sample, "/", sample, "-nanostringCNV-T.RCC")
gfile <- paste0(prefix, "/", sample, "/", sample, "-nanostringCNV-B.RCC")

readRaw <- function(file) {
  dat <- read.csv2(file, skip=26, sep=',', stringsAsFactors=F)
  dat <- dat[1:(dim(dat)[1]-3),]

  endo <- subset(dat, CodeClass=='Endogenous')
  inv <- subset(dat, CodeClass=='Invariant')
  ctrls <- subset(dat, CodeClass=='Positive' | CodeClass=='Negative')
  res <- subset(dat, CodeClass=='RestrictionSite')

  return(list(endo=endo,
              inv=inv,
              ctrls=ctrls,
              res=res))
}

## NORMALIZATION
# +ve control
posCtrls <- sapply(list.files(path = prefix, pattern = "nanostringCNV.*\\.RCC$", recursive = TRUE, full.names = TRUE), function(i) subset(readRaw(i)$ctrls, CodeClass=="Positive")$Count)
geomSample <- apply(posCtrls, 2, mean)
geomAll <- exp(mean(log(geomSample)))
sf <- geomSample / geomAll

# -ve control
negCtrls <- sapply(list.files(path = prefix, pattern = "nanostringCNV.*\\.RCC$", recursive = TRUE, full.names = TRUE), function(i) subset(readRaw(i)$ctrls, CodeClass=="Negative")$Count)
meanSample <- apply(negCtrls, 2, mean)
sdSample <- apply(negCtrls, 2, sd)
bg <- meanSample + 3 * sdSample

# calculating mean INV count across all samples
randINVmeancount <- mean(sapply(list.files(path = prefix, pattern = "nanostringCNV.*\\.RCC$", recursive = TRUE, full.names = TRUE), function(i) mean(readRaw(i)$inv$Count)))

# normalizing endogenous counts based on +ve control, -ve control and invariant counts
normNS <- function(NSdat, NSfile) {
  en <- NSdat$endo
  en$Countpos <- NSdat$endo$Count * sf[NSfile]
  en$Countneg <- en$Countpos - bg[NSfile]
  invnf <- randINVmeancount / mean(NSdat$inv$Count)
  en$Count <- en$Countneg* invnf
  return(en)
}

normNSinvar <- function(NSdat, NSfile) {
  invar <- NSdat$inv
  invar$Countpos <- NSdat$inv$Count * sf[NSfile]
  invar$Countneg <- invar$Countpos - bg[NSfile]
  invnf <- randINVmeancount / mean(NSdat$inv$Count)
  invar$Count <- invar$Countneg* invnf
  return(invar)
}

# normalization across probes
averageNS <- function(NSdat) {
  return(ddply(NSdat$endo_norm, ~Accession,summarize,mean=mean(Count)))
}

# tumor/normal ratio per probe
ratioNS <- function(tumor, normal) {
  ratio <- (tumor$average$mean / normal$average$mean) * 2
  names(ratio) <- tumor$average$Accession
  return(ratio)
}

# read in raw RCC files
tdat <- readRaw(tfile)
gdat <- readRaw(gfile)

# invariant normalization
tdat$endo_norm <- normNS(tdat, tfile)
gdat$endo_norm <- normNS(gdat, gfile)

# normalization across probes
tdat$average <- averageNS(tdat)
gdat$average <- averageNS(gdat)

final <- ratioNS(tdat, gdat)
final.df <- data.frame(Gene=names(final),Count=as.numeric(final))

# round (0.0 - 0.4) to 0 and (0.6 - 1.0) to 1 for Copy Number estimates
final.df$CopyNumber <- ifelse(signif(final.df$Count,2)%%0.5 == 0, signif(final.df$Count,2), round_any(final.df$Count,1))

write.table(final.df, file=paste0(prefix, "/", sample, "/", sample, ".nanostringCNV.out"), sep="\t", row.names=FALSE, quote=FALSE)
