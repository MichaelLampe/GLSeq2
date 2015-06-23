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
# Prints out a check to the log file if the counting file is correct.
#
header.message <- paste("echo \"Checking SAM files for RSEM compatability. The file is: ",countable.sam,"\"",sep="")
header.file <- paste(header.message,">> ",log.file)
check.file <- paste("rsem-sam-validator",countable.sam,">>",log.file)
#
refName <- paste(this.resName,"RSEM","index",sep=".")
prepareReference <- paste("rsem-prepare-reference",refFASTAname,refName)
#
countFile <- paste(this.resName,"RSEM","counts",sep=".")
rsemOptions <- paste("--calc-pme","--calc-ci","-p",nCores)
if (paired.end) rsemOptions <- paste(rsemOptions,"--sam","--paired-end")
if (!paired.end) rsemOptions <- paste(rsemOptions,"--sam")
#
if (!is.na(countable.sam)){
  calculateExpression <- paste("rsem-calculate-expression",rsemOptions,countable.sam,refName,countFile)
  # Movement
  if (count.comm != "") count.comm <- paste(count.comm,"&&",check.file,"&&", prepareReference, "&&",calculateExpression)
  if (count.comm == "") count.comm <- paste(header.file,"&&",check.file,"&&",prepareReference, "&&", calculateExpression)
}