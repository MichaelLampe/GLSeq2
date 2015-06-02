source ("GLSeq.Util.R")
#
#############
# Null Fixes
#############
#
if (is.null(this.resName)) this.resName <- text.add
#
#
refName <- paste(this.resName,"index",sep=".")
prepareReference <- paste("rsem-prepare-reference",refFASTAname,refName)

countfile <- paste(this.resName,"RSEM","counts", sep=".") 
rsemOptions <- paste("--sam","--paired-end")
calculateExpression <- paste("rsem-calculate-expression",rsemOptions,countable.sam,refName,countFile)
if (count.comm != "") count.comm <- paste(count.comm, "&&", prepareReference, "&&",calculateExpression)
if (count.comm == "") count.comm <- paste(prepareReference, "&&", calculateExpression)

