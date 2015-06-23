args <- commandArgs(trailingOnly = TRUE)
countable.sam <- as.character(args[1])
rGenome <- as.character(args[2])
refGFFname <- as.character(args[3])
dest.dir <- as.character(args[4])
this.resName <- as.character(args[5])
paired.end <- as.logical(args[6])

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
if (paired.end)try(featureCounts(files=countable.sam,annot.ext=refGFFname,isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id",isPairedEnd=TRUE,nthreads=1))
if(!paired.end)try(featureCounts(files=countable.sam,annot.ext=refGFFname,isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id",isPairedEnd=FALSE,nthreads=1))
sink()