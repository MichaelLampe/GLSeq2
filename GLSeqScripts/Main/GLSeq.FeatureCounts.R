args <- commandArgs(trailingOnly = TRUE)
countable.sam <- as.character(args[1])
rGenome <- as.character(args[2])
refGFFname <- as.character(args[3])
dest.dir <- as.character(args[4])
this.resName <- as.character(args[5])
paired.end <- as.logical(args[6])
idAttr <- as.character(args[7])

library(Rsubread)
#############
# Makes sure in the correct dir
#############
setwd(dest.dir)
#
#
#
#############
# Reads output into a file
#############
countfile <- paste(this.resName,"FeatureCounts","counts", sep=".")
sink(file = countfile,append=TRUE,type=c("output"),split=FALSE)
featureCounts(files=countable.sam,annot.ext=paste(dest.dir,refGFFname),isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType=idAttr,isPairedEnd=paired.end,nthreads=1)
sink()