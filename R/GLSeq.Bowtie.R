source("GLSeq.Util.R")
setwd(dest.dir)
# copy genome indices to the destimation dir: 
ref.dir <- paste(base.dir, rGenome, sep="")
indCopy <- paste("cd ", ref.dir, " && cp * ", dest.dir, sep="")
system(indCopy)
#
comm.stack.pool <- NULL # 
#
# Argument the form:
# Prepare reference, what aligner to use, what reference to use, file name
#
#-f means FASTA input, --seed indicates random number generator seed will stay constant.
IndexOptions <- paste("-f","--seed 0")
# --sam gives a sam file output, -t gives time information for the process.
BowtieOptions <- paste("-S","-t","-q")
Bowtie2Options <- paste ("-t")
if (qAlgor == "Bowtie"){
  index <- paste("bowtie-build",IndexOptions,refFASTAname,rGenome)
}
if (qAlgor == "Bowtie2"){
  index <- paste("bowtie2-build",IndexOptions,refFASTAname,rGenome)
}
comm.stack.pool <- paste(index)
for (zz in 1:nStreams) {
  # assembly and runing the system command, one library at a time:
  for (i in rangelist[[zz]]) {
    ###################
    # Alignment with SAM output
    ###################
    # names of current fastq files:
    fq.left <- fqfiles.table[i,1]
    if (paired.end) fq.right <- fqfiles.table[i,2]
    this.library <- substr(fqfiles.table[i,1], 1,libNchar)
    this.resName <- paste(this.library, text.add, sep=".")
    unsorted.sam <- paste(this.resName, "unsorted.sam", sep=".")
    ###################
    # Alignment
    ###################
    if (qAlgor == "Bowtie"){
      if (paired.end){
        align <- paste("bowtie",BowtieOptions,rGenome,"-1",fq.left,"-2",fq.right,unsorted.sam)
      }else{
        align <- paste("bowtie",BowtieOptions,rGenome,fq.left,unsorted.sam)
      }
    }
    if (qAlgor == "Bowtie2"){
      if (paired.end){
        align <- paste("bowtie2","-fr",Bowtie2Options,rGenome,"-U","-1",fq.left,"-2",fq.right,"-S",unsorted.sam)
      }else{
        align <- paste("bowtie2",Bowtie2Options,rGenome,"-U",fq.left,"-S",unsorted.sam)
      }
    }
    ###################
    # Counting Step
    ###################
    #
    # TODO: Need to add cleanup steps related to counting
    #
    if (!is.null(cAlgor)){
      source("GLSeq.Counting.R")
    }
    if (!is.null(count.comm)){
      comm.stack.pool <- paste(comm.stack.pool " && ",align," && ", count.comm,sep="")
    } else{
    comm.stack.pool <- paste(comm.stack.pool " && ",align,sep="")
    }
  }
  if (zz ==1) fileCompletenessID <- paste(text.add, ".completeExpression", sep="")
  comm.stack.pool <- paste(comm.stack.pool,  " && echo  >  ", fileCompletenessID, ".", zz, " & ", sep="")
}
# collection of the results: 
collLog <- paste(destDirLog, text.add, ".ResultsCollectLog.txt", sep="")
collerr <- paste(destDirLog, text.add, ".ResultsCollectErrors.txt", sep="")
collResults <- paste("cd ", base.dir, " && ", "Rscript GLSeqResultsCollect.R ", text.add, base.dir, dest.dir, " 1>> ", collLog, " 2>> ", collerr, " &", sep="") 
if (resCollect == "nocollect") collResults <- "\n"
# 
# pool of all the system commands (versions for +/- compute expression and +/- collect results): 
#
comm.stack.pool <- paste(comm.stack.pool, collResults, "\n", sep=" ") 
if (exprRun == "noexprcalc") comm.stack.pool <- paste(collResults, "\n", sep=" ") 
  