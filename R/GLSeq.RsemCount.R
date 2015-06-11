source ("GLSeq.Util.R")
#
#
setwd(destDirRSEMCount)
#############
# Null Fixes
#############
#
if (is.null(this.resName)) this.resName <- text.add
#
#
refName <- paste(this.resName,"index",sep=".")
prepareReference <- paste("rsem-prepare-reference",refFASTAname,refName)
#
countFile <- paste(this.resName,"RSEM","counts",sep=".")
if (paired.end) rsemOptions <- paste("--sam","--paired-end")
if (!paired.end) rsemOptions <- paste("--sam")
#
if (!is.na(countable.sam)){
  calculateExpression <- paste("rsem-calculate-expression",rsemOptions,countable.sam,refName,countFile)
  # Movement
  if (count.comm != "") count.comm <- paste(count.comm, "&&", prepareReference, "&&",calculateExpression)
  if (count.comm == "") count.comm <- paste(prepareReference, "&&", calculateExpression)
}