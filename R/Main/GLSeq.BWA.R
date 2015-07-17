source("GLSeq.Util.R")
source("GLSeq.Alignment.Functions.R")
setwd(dest.dir)

comm.stack.pool <- NULL

indCopy <- copyGenome(base.dir,rGenome,refFASTAname,dest.dir)
system(indCopy)
###################################################################################
# Index the BWA 
###################################################################################
index <- paste(bwaPath, "index", refFASTAname)
system(index)
#
for (zz in 1:nStreams) {
  # Assembly and runing the system command, one library at a time:
  for (i in rangelist[[zz]]) {
    ###################
    # Alignment with SAM output
    ###################
    # Grab and name everything correctly
    fq.left <- fqfiles.table[i,1]
    if (paired.end) fq.right <- fqfiles.table[i,2]
    name <- assign.name(fqfiles.table[i,1],paired.end)
    this.resName <- assign.resName(name,text.add)
    # Names of the expected sai files
    # Sai files are an intermediate file type for BWA only.
    # Read the original paper for more info: http://bioinformatics.oxfordjournals.org/content/25/14/1754.full
    sainame.left <- paste(fq.left,"sai",sep=".")
    if (paired.end) sainame.right <- paste(fq.right,"sai",sep=".")
    # Alignment Commands
    # Two alignments happen no matter what, after the SAI file is made we can combine paired files into one SAM
    # file using some BWA tools.
    aln.left <- paste(bwaPath, "aln", refFASTAname, fq.left, ">", sainame.left)
    if (paired.end) aln.right <- paste(bwaPath, "aln", refFASTAname, fq.right, ">", sainame.right)
    # Sam File Creation
    # Sampe/Samse transforms a .sai file into a SAM file.
    # We take the output of the Sam*e program and appends the content to the unsorted SAM file
    if (paired.end) sam.create <- paste(bwaPath, "sampe", refFASTAname, sainame.left, sainame.right, fq.left, fq.right, ">>", unsorted.sam)
    if (!(paired.end)) sam.create <- paste(bwaPath, "samse", refFASTAname, sainame.left, fq.left, ">>", unsorted.sam)
    #
    ###################
    # SAM file Cleanup
    ###################
    # Picard tools has a Sam cleaner that helps us out here
    cleanSAM <- paste("java -Xmx2g -jar ",picardToolsPath, "CleanSam.jar", sep="")
    # We give the file that will be output from the CleanSam.jar the Resname + a cleaned.sam suffix
    cleaned.sam <- paste(this.resName, "cleaned.sam", sep=".")
    # 
    # SAM cleanup system command:
    # I = Input file; O= Output file
    cleansam.comm <- paste(cleanSAM, " I=", unsorted.sam, " O=", cleaned.sam, sep="") # System command #4
    #
    ###################
    # Adding RG Header + sorting
    ###################
    # Picard tool also lets us modifty the Read group headers
    headersortSAM <- paste("java -Xmx2g -jar ",picardToolsPath, "AddOrReplaceReadGroups.jar", sep="")
    # Name of the processed (final) SAM file
    final.sam <- paste(this.resName, "final.sam", sep=".")
    # 
    # Implement the logic for sequencer specificity here.  
    #
    # This solution may not work for data derived from Non-illumina sequencing data. Haven't tested that yet as I don't have any data to.
    # I = Input file; O= Output File
    finalsam.comm <- paste(headersortSAM, " I=", cleaned.sam, " O=", final.sam, " SO=coordinate LB=", refFASTAname, " PL=ILLUMINA PU=unknown SM=", this.resName, " VALIDATION_STRINGENCY=LENIENT", sep="")
    ###################
    # SAM => BAM file with index 
    ###################
    # Just some name things to get us through the Samtools sort.
    sorted.arg <- paste(this.resName, "sorted", sep=".")
    ref.index <- paste(refFASTAname, "fai", sep=".") 
    sorted.bam <- paste(this.resName, "sorted.bam", sep=".") 
    # -u = Uncompressed BAM file output, better for the pipe
    # -S = Input is a SAM File
    # -t = TAB-delimited file
    # I'm not sure if we need to first pipe it from into the sorted.arg format and index or if we could just sort
    # This is something I will test in the future to possibly speed up the protocol a bit.
    bam.create <- paste("samtools view -uS -t ", ref.index, final.sam, " | samtools sort - ", sorted.arg)
    bam.index <- paste("samtools index", sorted.bam) # System command #7
    ###################
    # Converting the final bam (coordinate-sorted)
    # to name-sorted  sam:  
    ###################
    # Countable.Sam is conserved among all the counting protocols, so the name of it is given via a function
    # For consistency as we move towards a counting protocol
    countable.sam <- countable.sam.name(this.resName)
    paired.arg <- paste(this.resName, "paired", sep=".")
    paired.bam <- paste(this.resName, "paired.bam", sep=".")
    # let's avoid pipe here to save RAM
    # Directly write the samtools output to a new countable file after sorting
    # -n in the samtools sort indicates that the document will be sorted by read name rather than chromosome
    # -h in the samtools view indicates a header will be printed
    countable.comm <- paste("samtools sort -n", sorted.bam, paired.arg, "&&", "samtools view -h", paired.bam, ">", countable.sam)
    ###################
    # Counting Step
    ###################
    #
    count.comm <- ""
    if (counting == "counting"){
      setwd(base.dir)
      source("GLSeq.Counting.R")
    }
    ###################
    # Housekeeping 
    ###################
    # We remove a bunch of files here to make sure we don't clog up the hard drives too much with all the intermediate steps
    if (paired.end) spaceCleanup <- paste("rm", sainame.left, "&& rm", sainame.right, "&& rm", unsorted.sam, "&& rm", cleaned.sam, "&& rm", final.sam, "&& rm", paired.bam)
    if (!(paired.end)) spaceCleanup <- paste("rm", sainame.left, "&& rm", unsorted.sam, "&& rm", cleaned.sam, "&& rm", final.sam, "&& rm", paired.bam)
    # Paired Ended Samples vs Unpaired
    # Paired end align two files are once (Because they are paired)
    # So we add in both teh right and the left.
      if (paired.end){
        comm.i <- paste(aln.left,"&&",aln.right)
      }
      if (!paired.end){
        comm.i <- paste(aln.left)
      }
      comm.i <- paste(comm.i,"&&",sam.create,"&&",cleansam.comm,"&&",finalsam.comm,"&&",bam.create,"&&",bam.index,"&&",countable.comm)
    # Counting
    # Whatever was added to count.comm when it looked for counting protocols is added here
    comm.i <- paste(comm.i,"&&",count.comm)
    # Cleanup File 
    # Removes unneeded files.
    comm.i <- paste(comm.i,"&&",spaceCleanup)
    # For the very first assembly in the stack: 
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
  comm.stack.pool <- paste(comm.stack.pool,"&")
} 
comm.stack.pool <- paste(comm.stack.pool,"wait")