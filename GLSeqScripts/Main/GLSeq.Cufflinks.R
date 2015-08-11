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
#############
# Null Fixes
#############
#
if (is.null(this.resName)) this.resName <- paste(dest.dir,text.add,sep="")
if(Condor){
  # We've allocated 8 cores in the Condor wrapper
  nCores <- 8
}
#
if (!is.na(countable.sam)){
  individual.run.dir <- paste(fqfiles.table[i,1],"IndividualCountData",sep=".")
  cufflinks.dest <- paste(destDirCufflinksCount,individual.run.dir,sep="")
  cufflinks.sorted.name <- paste("cufflinks.ready",countable.sam,sep=".")
  #
  create.run.dir <- paste("mkdir",cufflinks.dest)
  cufflinks.unsorted.name <- paste(countable.sam,"cufflinks.unsorted.sam",sep=".")
  # Cufflinks requires coordinate sorted.
  cufflinks.ready.create <- paste("samtools view -uS",countable.sam,"-o",cufflinks.unsorted.name)
  if (Condor) {
    cufflinks.ready.sort <- paste("samtools sort -@ 6 -m 32G", cufflinks.unsorted.name, cufflinks.sorted.name)
  } else {
    cufflinks.ready.sort <- paste("samtools sort ", cufflinks.unsorted.name, cufflinks.sorted.name)
  }
  # This will end up being a BAM file due to how Samtools sorts things.
  cufflinks.ready.countable <- paste(cufflinks.sorted.name,"bam",sep=".")
  cufflinks.ready.countable <- paste(dest.dir,cufflinks.ready.countable,sep="")
  #
  # No progress bar makes everything feel better.
  cufflinks.options <- paste("--quiet","-p",nCores,"--output-dir",cufflinks.dest)
  cufflinks.run <- paste("cufflinks",cufflinks.options,cufflinks.ready.countable)
  #
  if (count.comm != "") count.comm <- paste(count.comm,"&&",create.run.dir,"&&",cufflinks.ready.create,"&&",cufflinks.ready.sort,"&&",cufflinks.run)
  if (count.comm == "") count.comm <- paste(create.run.dir,"&&",cufflinks.ready.create,"&&",cufflinks.ready.sort,"&&",cufflinks.run)
}
