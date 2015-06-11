#########################################################
# Great Lakes Seq package for low-level processing of RNA-Seq data
# Trimmomatic-BWA-HTSeq quantification of the expression values
# TopHat
# June 2015
#########################################################
#
#
#
source("GLSeq.Util.R")
setwd(dest.dir)
####################################
# Copy genome indices to the destimation dir: 
####################################
ref.dir <- paste(base.dir, rGenome, sep="")
indCopy <- paste("cd ", ref.dir, " && cp ",refFASTAname," ",dest.dir, sep="")
system(indCopy)
#
comm.stack.pool <- NULL # 

for (zz in 1:nStreams) {
  for (i in rangelist[[zz]]) {
    ###################
    # Alignment with SAM output
    ###################
    # names of current fastq files:
    fq.left <- fqfiles.table[i,1]
    if (paired.end) fq.right <- fqfiles.table[i,2]
    this.library <- substr(fqfiles.table[i,1], 1,libNchar)
    this.resName <- paste(this.library, text.add, sep=".")
    unsorted.sam <- paste(this.resName, "unsorted", sep=".")
    
    
    
  }
}