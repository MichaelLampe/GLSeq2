source("GLSeq.Util.R")
source("GLSeq.Alignment.Functions.R")
setwd(dest.dir)

indCopy <- copy.genome(base.dir,rGenome,refFASTAname,dest.dir)
system(indCopy)
###################################################################################
# Index the BWA 
###################################################################################
index <- paste(bwaPath, "index", refFASTAname)
system(index)
#
for (zz in 1:nStreams) {
  # Assembly and runing the system command, one library at a time:
  # For the very first assembly in the stack: 
  if (zz == 1) comm.stack.pool <- "date"
  if (zz != 1) comm.stack.pool <- paste(comm.stack.pool,"date")
  
  for (i in rangelist[[zz]]) {
    # Names of current fastq files:
    fq.left <- fqfiles.table[i,1]
    if (paired.end) fq.right <- fqfiles.table[i,2]
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
    aln.left <- paste(bwaPath, "aln", refFASTAname, fq.left, ">", sainame.left)
    if (paired.end) aln.right <- paste(bwaPath, "aln", refFASTAname, fq.right, ">", sainame.right)
    # Sam File Creation
    # Sampe/Samse transforms a .sai file into a SAM file.
    # We take the output of the Sam*e program and appends the content to the unsorted SAM file
    if (paired.end) sam.create <- paste(bwaPath, "sampe", refFASTAname, sainame.left, sainame.right, fq.left, fq.right, ">>", unsorted.sam)
    if (!(paired.end)) sam.create <- paste(bwaPath, "samse", refFASTAname, sainame.left, fq.left, ">>", unsorted.sam)
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
    bam.create <- paste("samtools view -uS -t ", ref.index, unsorted.sam, " | samtools sort - ", sorted.arg)
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
    bai <- paste(sorted.bam,"bai",sep=".")
    if (paired.end) spaceCleanup <- paste("rm",fq.left,"&& rm",fq.right,"&& rm", sainame.left, "&& rm", sainame.right, "&& rm", unsorted.sam, "&& rm", paired.bam,"&& rm",sorted.bam,"&& rm", bai)
    if (!(paired.end)) spaceCleanup <- paste("rm", sainame.left, "&& rm",fq.left,"&& rm", unsorted.sam,"&& rm", paired.bam,"&& rm",sorted.bam,"&& rm", bai)
    # Paired Ended Samples vs Unpaired
    # Paired end align two files are once (Because they are paired)
    # So we add in both teh right and the left.
      if (paired.end){
        comm.i <- paste(aln.left,"&&",aln.right)
      }
      if (!paired.end){
        comm.i <- paste(aln.left)
      }
      comm.i <- paste(comm.i,"&&",sam.create,"&&",bam.create,"&&",bam.index,"&&",countable.comm)
    # Counting
    # Whatever was added to count.comm when it looked for counting protocols is added here
    comm.i <- paste(comm.i,"&&",count.comm)
    # Cleanup File 
    # Removes unneeded files.
    comm.i <- paste(comm.i,"&&",spaceCleanup)
    # system(comm.i)
    # For the very first assembly in the stack: 
    # For the very first assembly in the stack (i = 1)
    comm.stack.pool <- paste(comm.stack.pool,"&&",comm.i)
  } # for i 
  comm.stack.pool <- paste(comm.stack.pool,"&")
} 
comm.stack.pool <- paste(comm.stack.pool,"wait")