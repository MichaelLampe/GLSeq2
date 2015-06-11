source("GLSeq.Util.R")
setwd(dest.dir)
# if (!(paired.end)) stop('BWA pipeline is currently supported for paired-end libraries only \n')
# copy genome indices to the destimation dir: 
ref.dir <- paste(base.dir, rGenome, sep="")
indCopy <- paste("cd ", ref.dir, " && cp ",refFASTAname," ",dest.dir, sep="")
system(indCopy)
#
comm.stack.pool <- NULL # 
###################################################################################
# Index the BWA 
###################################################################################
# Might be able to do this with rsem too in some cases (Create index file)
###################################################################################
index <- paste(bwaPath, "index", refFASTAname) # System command #0
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
    unsorted.sam <- paste(this.resName, "unsorted.sam", sep=".")
    #
    # names of the expected sai files:
    sainame.left <- paste(fq.left,"sai",sep=".")
    if (paired.end) sainame.right <- paste(fq.right,"sai",sep=".")
    # Alignment Commands:
    aln.left <- paste(bwaPath, "aln", refFASTAname, fq.left, ">", sainame.left) # System command #1
    if (paired.end) aln.right <- paste(bwaPath, "aln", refFASTAname, fq.right, ">", sainame.right) # System command #2
    # Sam File Creation:
    if (paired.end) sam.create <- paste(bwaPath, "sampe", refFASTAname, sainame.left, sainame.right, fq.left, fq.right, ">>", unsorted.sam) # System command #3
    if (!(paired.end)) sam.create <- paste(bwaPath, "samse", refFASTAname, sainame.left, fq.left, ">>", unsorted.sam)
    #
    ###################
    # SAM file Cleanup
    ###################
    # the executable:
    cleanSAM <- paste("java -Xmx2g -jar ",picardToolsPath, "CleanSam.jar", sep="")
    # name of the cleaned SAM file:
    cleaned.sam <- paste(this.resName, "cleaned.sam", sep=".")
    # 
    # SAM cleanup system command:
    # I = Input file; O= Output file
    cleansam.comm <- paste(cleanSAM, " I=", unsorted.sam, " O=", cleaned.sam, sep="") # System command #4
    #
    ###################
    # Adding RG Header + sorting
    ###################
    # the executable: 
    headersortSAM <- paste("java -Xmx2g -jar ",picardToolsPath, "AddOrReplaceReadGroups.jar", sep="")
    # name of the processed (final) SAM file:
    final.sam <- paste(this.resName, "final.sam", sep=".")
    # 
    # Implement the logic for sequencer specificity here.  
    #
    finalsam.comm <- paste(headersortSAM, " I=", cleaned.sam, " O=", final.sam, " SO=coordinate LB=", refFASTAname, " PL=ILLUMINA PU=unknown SM=", this.resName, " VALIDATION_STRINGENCY=LENIENT", sep="") # System command #5
    #
    ###################
    # SAM => BAM file with index 
    ###################
    sorted.arg <- paste(this.resName, "sorted", sep=".")
    ref.index <- paste(refFASTAname, "fai", sep=".") 
    sorted.bam <- paste(this.resName, "sorted.bam", sep=".") 
    # -u = Uncompressed BAM file output, better for the pipe
    # -S = Input is a SAM File
    # -t = TAB-delimited file
    bam.create <- paste("samtools view -uS -t ", ref.index, final.sam, " | samtools sort - ", sorted.arg) # System command #6
    bam.index <- paste("samtools index", sorted.bam) # System command #7
    #
    ###################
    # Converting the final bam (coordinate-sorted)
    # to name-sorted  sam:  
    ###################
    #
    countable.sam <- paste(this.resName, "countable.sam", sep=".")
    paired.arg <- paste(this.resName, "paired", sep=".")
    paired.bam <- paste(this.resName, "paired.bam", sep=".")
    # let's avoid pipe here to save RAM
    # Directly write the samtools output to a new countable file after sorting
    # -n in the samtools sort indicates that the document will be sorted by read name rather than chromosome
    # -h in the samtools view indicates a header will be printed
    countable.comm <- paste("samtools sort -n", sorted.bam, paired.arg, "&&", "samtools view -h", paired.bam, ">", countable.sam) # System command #8
    ###################
    # Counting Step
    ###################
    #
    count <- FALSE
    count.comm <- paste ("")
    if (counting == "counting"){
      setwd(base.dir)
      count <- TRUE
      source("GLSeq.Counting.R")
    }
    ###################
    # Housekeeping 
    ###################
    if (paired.end) spaceCleanup <- paste("rm", sainame.left, "&& rm", sainame.right, "&& rm", unsorted.sam, "&& rm", cleaned.sam, "&& rm", final.sam, "&& rm", paired.bam) # System command #10
    if (!(paired.end)) spaceCleanup <- paste("rm", sainame.left, "&& rm", unsorted.sam, "&& rm", cleaned.sam, "&& rm", final.sam, "&& rm", paired.bam) # System command #10
    #
    # Paired Ended Samples vs Unpaired
    #
      if (paired.end){
        comm.i <- paste(aln.left,"&&",aln.right)
      }
      if (!paired.end){
        comm.i <- paste(aln.left)
      }
      comm.i <- paste(comm.i,"&&",sam.create,"&&",cleansam.comm,"&&",finalsam.comm,"&&",bam.create,"&&",bam.index,"&&",countable.comm)
    #
    # Counting
    # 
    if (count) comm.i <- paste(comm.i,"&&",count.comm)
    #
    # Cleanup File 
    #
    comm.i <- paste(comm.i,"&&",spaceCleanup)
    #
    # for the very first assembly in the stack: 
    if (i == rangelist[[zz]][1])  comm.stack.pool <- paste(comm.stack.pool, "date && ", comm.i)
    # for subsequent assemblies of every stack: 
    if (i != rangelist[[zz]][1])  comm.stack.pool <- paste(comm.stack.pool, " && date && ", comm.i)
    # system(comm.i)
    if (resCollect == "collect"){
      collLog <- paste(destDirLog, text.add, ".ResultsCollectLog.txt", sep="")
      collerr <- paste(destDirLog, text.add, ".ResultsCollectErrors.txt", sep="")
      collResults <- paste("cd ", base.dir, " && ", "Rscript GLSeqResultsCollect.R ", text.add, base.dir, dest.dir, " 0 1>> ", collLog, " 2>> ", collerr, sep="")
      if (is.null(comm.stack.pool)) comm.stack.pool <- paste(collResults)
      if (!is.null(comm.stack.pool)) comm.stack.pool <- paste(comm.stack.pool,"&&",collResults)
    }
  } # for i 
  if (zz ==1) fileCompletenessID <- paste(text.add, ".completeExpression", sep="")
  comm.stack.pool <- paste(comm.stack.pool,  " && echo  >  ", fileCompletenessID, ".", zz, " & ", sep="")
} # for zz