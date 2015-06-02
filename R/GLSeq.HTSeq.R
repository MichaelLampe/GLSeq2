source("GLSeq.Util.R")
if (is.null(countable.sam)){
  countable.sam <- paste(this.resName, "countable.sam", sep=".")
}
countfile <- paste(this.resName, "counts", sep=".") 
# initializing strandness option with a non-strand-specific value:
countOpt <-  "--stranded=no"
# overwriting countOpt with information on strand-specificity:
if (libstrand == "R") countOpt <-  "--stranded=reverse"
if (libstrand == "F") countOpt <-  "--stranded=yes"
# selecting ID attribute for counts reporting: 
countOpt <- paste(countOpt, " --idattr=", idAttr, sep="")
count.comm <- paste("python -m HTSeq.scripts.count", countOpt, countable.sam, refGFFname, " > ", countfile)