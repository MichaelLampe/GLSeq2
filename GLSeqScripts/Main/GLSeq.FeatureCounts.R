args <- commandArgs(trailingOnly = TRUE)
countable.sam <- as.character(args[1])
destDirFeatureCounts <- as.character(args[2])
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
#############
# Reads output into a file
#############
counts.summary <- paste(this.resName,".FeatureCounts.summary.txt",sep="")
counts.data <- paste(this.resName,".FeatureCounts.counts.csv",sep="")
counts.stats <- paste(this.resName,".FeatureCounts.stats.csv",sep="")
counts.annotation <- paste(this.resName,".FeatureCounts.annotations.csv",sep="")
sink(file = counts.summary,append=TRUE,type=c("output"),split=FALSE)
fc <- featureCounts(files=countable.sam,annot.ext=paste(dest.dir,refGFFname,sep=""),isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id",isPairedEnd=paired.end,nthreads=1)
sink()
write.csv(fc[1],counts.data)
write.csv(fc[2],counts.annotation)
write.csv(fc[4],counts.stats)