source("GLSeq.Util.R")
setwd(dest.dir)
# copy genome indices to the destimation dir: 
ref.dir <- paste(base.dir, rGenome, sep="")
indCopy <- paste("cd ", ref.dir, " && cp ",refFASTAname," ",dest.dir, sep="")
system(indCopy)
#
comm.stack.pool <- NULL # 
####################
# Index the Bowtie or Bowtie2 Aligner
#
#Argument the form:
# Prepare reference, what aligner to use, what reference to use, file name
#
#-f means FASTA input
####################
IndexOptions <- paste("-f")
if (aAlgor == "Bowtie"){
  index <- paste("bowtie-build",IndexOptions,refFASTAname,rGenome)
}
if (aAlgor == "Bowtie2"){
  index <- paste("bowtie2-build",IndexOptions,refFASTAname,rGenome)
}
system(index)
#
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
    countable.sam <- paste(this.resName, "countable.sam", sep=".")
    ###################
    # Alignment
    ###################
    # --sam gives a sam file output, -t gives time information for the process.
    # --chunkmbs may fix paired end read problems skipping reads due to memory exhaustion
    # Default for chunkmbs = 64, 512 is 8x larger.  No known problems created by this option.
    ###################
    if (aAlgor == "Bowtie"){
      BowtieOptions <- paste("-S","-t","-q","--chunkmbs","512")
      if (paired.end){
        align <- paste("bowtie",BowtieOptions,rGenome,"-1",fq.left,"-2",fq.right,countable.sam)
      }else{
        align <- paste("bowtie",BowtieOptions,rGenome,fq.left,countable.sam)
      }
    }
    ###################
    # Might need to change a few things in Bowtie2
    # It gets hung soemtimes, though I'm not sure 
    # if it is the same issue as in the Bowtie algo
    ###################
    if (aAlgor == "Bowtie2"){
      #
      # All these options are necessary just to get it to conform with the RSEM counting protocol... (Except -t and -q)
      # 
      Bowtie2Options <- paste("-t","-q","--sensitive","--dpad 0","--gbar 99999999","--mp 1,1","--score-min L,0,-0.1")
      if (paired.end){
        Bowtie2Options <- paste(Bowtie2Options,"--no-mixed","--no-discordant")
        align <- paste("bowtie2",Bowtie2Options,rGenome,"-1",fq.left,"-2",fq.right,"-S",countable.sam)
      }else{
        align <- paste("bowtie2",Bowtie2Options,rGenome,"-U",fq.left,"-S",countable.sam)
      }
    }
    ###################
    # Counting Step
    ###################
    #
    count <- FALSE
    if (counting == "counting"){
      setwd(base.dir)
      count <- TRUE
      source("GLSeq.Counting.R")
    }
    ###################
    # Command Construction
    ###################
    comm.i <- paste(align)
    if (count) comm.i <- paste(comm.i, "&&", count.comm)
    # for the very first assembly in the stack: 
    if (i == rangelist[[zz]][1])  comm.stack.pool <- paste(comm.stack.pool, " date && ", comm.i)
    # for subsequent assemblies of every stack: 
    if (i != rangelist[[zz]][1])  comm.stack.pool <- paste(comm.stack.pool, " && date && ", comm.i)
    #
    if (resCollect == "collect"){
      collLog <- paste(destDirLog, text.add, ".ResultsCollectLog.txt", sep="")
      collerr <- paste(destDirLog, text.add, ".ResultsCollectErrors.txt", sep="")
      collResults <- paste("cd ", base.dir, " && ", "Rscript GLSeqResultsCollect.R ", text.add, base.dir, dest.dir, " 0 1>> ", collLog, " 2>> ", collerr, sep="")
      if (is.null(comm.stack.pool)) comm.stack.pool <- paste(collResults)
      if (!is.null(comm.stack.pool)) comm.stack.pool <- paste(comm.stack.pool,"&&",collResults)
    }
    
  }
  if (zz ==1) fileCompletenessID <- paste(text.add, ".completeExpression", sep="")
  comm.stack.pool <- paste(comm.stack.pool,  " && echo  >  ", fileCompletenessID, ".", zz, " & ",sep="")
} 