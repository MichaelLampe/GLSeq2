#########################################################
# Great Lakes Seq package for low-level processing of RNA-Seq data
# Trimmomatic-BWA-HTSeq quantification of the expression values
# Cufflinks Counting Method
#
# Date: June 2015
# Author: Michael Lampe
#########################################################

source ("GLSeq.Util.R")
#
ref.dir <- paste(base.dir, rGenome, sep="")
refCopy <- paste("cd ", ref.dir, " && cp ",refGFFname," ",dest.dir, sep="")
system(refCopy)
setwd(dest.dir)
#############
# Null Fixes
#############
#
if (is.null(this.resName)) this.resName <- text.add
#
#
if (!is.na(countable.sam)){
  individual.run.dir <- paste(countable.sam,"IndividualCountData",sep=".")
  cufflinks.dest <- paste(destDirCufflinksCount,individual.run.dir,sep="")
  cufflinks.sorted.name <- paste("cufflinks.ready",countable.sam,sep=".")
  #
  create.run.dir <- paste("mkdir",cufflinks.dest)
  # Cufflinks requires coordinate sorted.
  cufflinks.ready.create <- paste("samtools view -uS",countable.sam, " | samtools sort - ", cufflinks.sorted.name)
  # This will end up being a BAM file due to how Samtools sorts things.
  cufflinks.ready.countable <- paste(cufflinks.sorted.name,"bam",sep=".")
  cufflinks.ready.countable <- paste(dest.dir,cufflinks.ready.countable,sep="")
  #
  # No progress bar makes everything feel better.
  cufflinks.options <- paste("--quiet --output-dir",cufflinks.dest)
  cufflinks.run <- paste("cufflinks",cufflinks.options,cufflinks.ready.countable)
  #
  if (count.comm != "") count.comm <- paste(count.comm,"&&",create.run.dir,"&&",cufflinks.ready.create,"&&",cufflinks.run)
  if (count.comm == "") count.comm <- paste(create.run.dir,"&&",cufflinks.ready.create,"&&",cufflinks.run)
}