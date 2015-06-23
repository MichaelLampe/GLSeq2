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
  cufflinks.sorted.name <- paste("cufflinks.ready",countable.sam,sep=".")
  cufflinks.ready.create <- paste("samtools view -uS",countable.sam, " | samtools sort - ", cufflinks.sorted.name)
  # This will end up being a BAM file due to how Samtools sorts things.
  cufflinks.ready.countable <- paste(cufflinks.sorted.name,"bam",sep=".")
  # Prints out a check to the log file if the counting file is correct.
  #
  cufflinks.options <- paste("--output-dir",destDirCufflinksCount)
  cufflinks.run <- paste("cufflinks",cufflinks.options,cufflinks.ready.countable)
  #
  if (count.comm != "") count.comm <- paste(count.comm,"&&",cufflinks.ready.countable,"&&",cufflinks.run)
  if (count.comm == "") count.comm <- paste(cufflinks.ready.countable,"&&",cufflinks.run)
}