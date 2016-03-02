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

# calculating mean INV count across all samples (change 'size' flag for random subset)
#randINVmeancount <- mean(sapply(sample(list.files(path = '/data/storage/patients/cnv/nanostring/platform_comp'), size = 94), function(i) mean(readRaw(i)$inv$Count)))
randINVmeancount <- 1350.169

# normalizing endogenous counts based on invariant counts
normNS <- function(NSdat) {
  nf <- randINVmeancount / mean(NSdat$inv$Count)
  en <- NSdat$endo
  en$Count <- NSdat$endo$Count * nf
  #invar <- NSdat$inv
  #invar$Count <- NSdat$inv$Count * nf
  return(en)
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
tdat$endo_norm <- normNS(tdat)
gdat$endo_norm <- normNS(gdat)

# normalization across probes
tdat$average <- averageNS(tdat)
gdat$average <- averageNS(gdat)

final <- ratioNS(tdat, gdat)
final.df <- data.frame(Gene=names(final),Count=as.numeric(final))
write.table(final.df, file=paste0(prefix, "/", sample, "/", sample, ".nanostringCNV.out"), sep="\t", row.names=FALSE, quote=FALSE)
