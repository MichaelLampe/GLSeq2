#########################################################
# Great Lakes Seq package for low-level processing of RNA-Seq data
# Trimmomatic-BWA-HTSeq quantification of the expression values
# HTSeq Counting Method
#
# Author: Michael Lampe
#########################################################

source("GLSeq.Util.R")
# Reference genome is necessary for counting
ref.dir <- paste(base.dir, rGenome, sep="")
refCopy <- paste("cd ", ref.dir, " && cp ",refGFFname," ",dest.dir, sep="")
system(refCopy)
#
setwd(dest.dir)
#
#############
# Null Fixes
#############
#
if (is.null(this.resName)) this.resName <- text.add
#
countfile <- paste(this.resName,"HTSeq","counts",sep=".") 
# initializing strandness option with a non-strand-specific value:
countOpt <-  "--stranded=no"
# overwriting countOpt with information on strand-specificity:
if (libstrand == "R") countOpt <-  "--stranded=reverse"
if (libstrand == "F") countOpt <-  "--stranded=yes"
# selecting ID attribute for counts reporting: 
countOpt <- paste(countOpt, " --idattr=", idAttr, sep="")
#
# Throws the output into a counting file and then moves the output over into the HTSeq folder.
HtSeq.comm <- paste("python -m HTSeq.scripts.count",countOpt, countable.sam, refGFFname, ">", countfile,"&&","mv",countfile,destDirHTSeqCount)
# Adds onto the end of or creates the count.comm which communicates the counting routine to the main command pool
if (count.comm != "") count.comm <- paste(count.comm, ";", HtSeq.comm)
if (count.comm == "") count.comm <- paste(HtSeq.comm)