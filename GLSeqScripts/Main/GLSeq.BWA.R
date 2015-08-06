source("GLSeq.Util.R")
source("GLSeq.Alignment.Functions.R")

indCopy <- copy.genome(base.dir,rGenome,refFASTAname,dest.dir)
printOrExecute(indCopy,Condor)
###################################################################################
# Index the BWA
###################################################################################
index <- paste(bwaPath, "index", paste(dest.dir,refFASTAname,sep="")
printOrExecute(index,Condor)
#
# Declare both pools as nulls
comm.stack.pool <- NULL
comm.stack.pools <- NULL
for (zz in 1:nStreams) {
  for (i in rangelist[[zz]]) {
    # Names of current fastq files:
    fq.left <- paste(dest.dir,fqfiles.table[i,1],sep="")
    if (paired.end) fq.right <- paste(dest.dir,fqfiles.table[i,2],sep="")
    name <- fq.left
    if (paired.end){
      name <- substr(name,1,nchar(name) - 5)
    } else{
      name <- substr(name,1,nchar(name) - 3)
    }
    this.resName <- assign.resName(name,text.add)
    unsorted.sam <- paste(this.resName, "unsorted", sep=".")
    # Names of the expected sai files
    # Sai files are an intermediate file type for BWA only.
    # Read the original paper for more info: http://bioinformatics.oxfordjournals.org/content/25/14/1754.full
    sainame.left <- paste(fq.left,"sai",sep=".")
    if (paired.end) sainame.right <- paste(fq.right,"sai",sep=".")
    # Alignment Commands
    # Two alignments happen no matter what, after the SAI file is made we can combine paired files into one SAM
    # file using some BWA tools.
    aln.left <- paste(bwaPath, "aln", paste(dest.dir,refFASTAname,sep=""), fq.left, ">", sainame.left)
    if (paired.end) aln.right <- paste(bwaPath, "aln", paste(dest.dir,refFASTAname,sep=""), fq.right, ">", sainame.right)
    # Sam File Creation
    # Sampe/Samse transforms a .sai file into a SAM file.
    # We take the output of the Sam*e program and appends the content to the unsorted SAM file
    if (paired.end) sam.create <- paste(bwaPath, "sampe", paste(dest.dir,refFASTAname,sep=""), sainame.left, sainame.right, fq.left, fq.right, ">>", unsorted.sam)
    if (!(paired.end)) sam.create <- paste(bwaPath, "samse", paste(dest.dir,refFASTAname,sep=""), sainame.left, fq.left, ">>", unsorted.sam)
    ###################
    # SAM => BAM file with index
    ###################
    # Just some name things to get us through the Samtools sort.
    sorted.arg <- paste(this.resName, "sorted", sep=".")
    unsorted.bam <- paste(this.resName,"unsorted.bam")
    # -u = Uncompressed BAM file output, better for the pipe
    # -S = Input is a SAM File
    # -t = TAB-delimited file
    # I'm not sure if we need to first pipe it from into the sorted.arg format and index or if we could just sort
    # This is something I will test in the future to possibly speed up the protocol a bit.
    bam.create <- paste("samtools view -uS", unsorted.sam,unsorted.bam)
    ###################
    # Converting the final bam (coordinate-sorted)
    # to name-sorted  sam:
    ###################
    # Countable.Sam is conserved among all the counting protocols, so the name of it is given via a function
    # For consistency as we move towards a counting protocol
    countable.sam <- countable.sam.name(this.resName)
    sorted.bam <- paste(this.resName, "sorted.bam", sep=".")
    # let's avoid pipe here to save RAM
    # Directly write the samtools output to a new countable file after sorting
    # -n in the samtools sort indicates that the document will be sorted by read name rather than chromosome
    # -h in the samtools view indicates a header will be printed
    sort.bam <- paste("samtools sort -n", unsorted.bam, sorted.bam)
    index.bam <- paste("samtools index", sorted.bam)
    countable.comm <- paste("samtools view -h", sorted.bam, ">", countable.sam)
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
    bai <- paste(sorted.bam,"bai",sep=".")
    # Paired Ended Samples vs Unpaired
    # Paired end align two files are once (Because they are paired)
    # So we add in both teh right and the left.
      if (paired.end){
        comm.i <- paste(aln.left,"&&",aln.right)
      }
      if (!paired.end){
        comm.i <- paste(aln.left)
      }
      comm.i <- paste(comm.i,"&&",sam.create,"&&",bam.create,"&&",sort.bam,"&&",index.bam"&&",countable.comm)
    # Counting
    # Whatever was added to count.comm when it looked for counting protocols is added here
    comm.i <- paste(comm.i,"&&",count.comm)
    #
    if (is.null(comm.stack.pool))
    comm.stack.pool <- paste(comm.stack.pool,"&&",comm.i)
  } # for i
  comm.stack.pool <- paste(comm.stack.pool,"&")
}
comm.stack.pool <- paste(comm.stack.pool,"wait")