#########################################################
# Great Lakes Seq package for low-level processing of RNA-Seq data
# Trimmomatic-BWA-HTSeq quantification of the expression values
# (paired end libraries); based on "Dana's script" + Irene input + Oleg's observations on the data
# (results of the latter: default headcrop is increased to 12 and Illumina artifacts removal is introduced)
# September 2013
#########################################################
#
#
source("GLSeq.Util.R")
setwd(dest.dir)
####################################
# Copy genome indices to the destimation dir: 
####################################
ref.dir <- paste(base.dir, rGenome, sep="")
indCopy <- paste("cd ", ref.dir, " && cp ",refFASTAname," ",dest.dir, sep="")
system(indCopy)
#
comm.stack.pool <- NULL # 
# This will index the CUSHAW with your FASTA file, necessary.  
# Previously the user needed to do this, but it should be easier if we just take care of
# This stuff on our end during each run (It is fast && cheap anyway)
index <- paste(CUSHAW.index.path,"-a","bwtsw",refFASTAname)
####################################
# ~Hot~ system command to index aligner
####################################
system(index)
####################################
### Let's speed this up.
### If using the GPU, alignment occurs in sequence so as to not overload the GPU memory (Each process is ~6GB)
### Our NVidia Titan can sometimes run two at a time without alignment corruption, but this is a safer route.
####################################
sam.create <- NULL
if (GPU.accel){
  for (zz in 1:nStreams) {
    for (i in rangelist[[zz]]) {
      ###################
      # Alignment with SAM output
      ###################
      # names of current fastq files:
      fq.left <- fqfiles.table[i,1]
      if (paired.end) fq.right <- fqfiles.table[i,2]
      this.library <- substr(fqfiles.table[i,1], 1,libNchar)
      this.resName <- paste(this.library, text.add, sep=".")
      unsorted.sam <- paste(this.resName, "unsorted", sep=".")
      #
      ###################
      # CUSHAW alignment:
      ###################
      ############
      # I've tried to increase the cores to nStreams*nCores 
      # (That should be how many cores the user has available), 
      # but it seems to just slow down process.  
      # With our Titan GPU, 6 cores is faster than 12 cores
      # However, running just the GPU a slower overall 
      # speed (No GPU speed dif + loss of work done by cores).  I've tested 5 and 7 
      # cores and 6 is the fastest between them. Thus, 6 appears to be a good max
      # number of cores.
      #
      # This result may not be globaly true for different
      # GPUs, but there is no plan to detect and react to the user's
      # GPU so this is the safest route to ensure reasonable speeds.
      # Importantly, the larger the files that are being aligned the
      # faster the GPU is compared to a single core (I've seen up to 50x)
      # Thus, the important case, where the calculation is difficult and lengthy
      # is best taken care of by limiting the number of cores the
      # user can use.
      # ~ Michael
      #
      ############
      if (nCores > 6) {
        nCores <- 6
        warning("Too many cores can slow alignment.  Adjusting maximum cores to 6.")
      }
      if (paired.end)  create <- paste(CUSHAW.GPU.path, "-r", refFASTAname, "-q", fq.right, fq.left, "-o", unsorted.sam,"-t", nCores)
      if (!(paired.end)) create <- paste(CUSHAW.GPU.path, "-r", refFASTAname, "-s", fq.left, "-o", unsorted.sam,"-t", nCores)
      # Checks to make sure that this process actually is GPU
      # The is.null is there for later on implementation of parallel GPU runs. (Hopefully!)
      if (is.null(sam.create)){
        sam.create <- paste(create)
      } else{
        sam.create <- paste(sam.create,"&&",create)
      }
      #
      # Fixes the order of processing so that
      # The SAM file gets created when the user
      # Preps data, but it still executes in sequence
      #
      GPUspecialCase <- TRUE
    } # i
  } #nStreams (zz)
} #GPU.Accel
if (!is.null(sam.create)) {
  Cushawgpu.special.case <- paste(sam.create)
}
if (GPUspecialCase) {
  system(Cushawgpu.special.case)
  warning("The alignment step has now completed.")
}
#################
# Rest of the cleanup to get to counting, or CPU only CUSHAW
#################
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
    unsorted.sam <- paste(this.resName, "unsorted", sep=".")
    #
    ###################
    # If no GPU, we can run all of the above files in parallel
    ###################
    if (!GPU.accel){
      if (paired.end) sam.create<- paste(CUSHAW.path, "-r", refFASTAname, "-q", fq.left, fq.right, "-o", unsorted.sam , "-t", nCores,"-t", nCores) 
      if (!paired.end) sam.create <- paste(CUSHAW.path, "-r", refFASTAname, "-s", fq.left, "-o", unsorted.sam,"-t", nCores)
    }
    ###################
    # Convert BAM to sorted BAM file
    ################### 
    sorted.arg <- paste(this.resName, "unsorted", sep=".")
    ref.index <- paste(refFASTAname, "fai", sep=".") 
    unsorted.bam <- paste(this.resName, "unsorted.bam", sep=".") 
    bam.create <- paste("samtools view -uS -t ", ref.index, unsorted.sam, " | samtools sort - ", sorted.arg) # System command #6
    bam.index <- paste("samtools index", unsorted.bam) # System command #7
    ###################
    # Sort bam file by name
    ################### 
    sorted <- paste(this.resName,"sorted",sep=".")
    sort.bam <- paste("samtools sort -n",unsorted.bam,sorted)
    ###################
    # Convert back to SAM file
    ################### 
    countable.sam <- paste(this.resName,"countable","sam",sep=".")
    # This removes certain read errors that can cause counting programs to crash
    remove <- paste("'!/\t\\*\t/'")
    sorted.bam <- paste(this.resName,"sorted.bam",sep=".")
    convert.to.sam <- paste ("samtools view -h",sorted.bam,"|","awk",remove,">",countable.sam)
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
    # Current command:
    ###################
    if (!GPU.accel) comm.i <- paste(sam.create, "&&",bam.create,"&&",bam.index,"&&",sort.bam,"&&",convert.to.sam)
    if (GPU.accel) comm.i <- paste(bam.create,"&&",bam.index,"&&",sort.bam,"&&",convert.to.sam)
    if (count.comm != "") comm.i <- paste(comm.i, "&&", count.comm)
    #
    # For the very first assembly in the stack (i = 1)
    if (i == rangelist[[zz]][1])  comm.stack.pool <- paste(comm.stack.pool, " date && ", comm.i)
    #
    # For subsequent assemblies of every stack (i > 1)
    if (i != rangelist[[zz]][1])  comm.stack.pool <- paste(comm.stack.pool, " && date && ", comm.i)
    # system(comm.i)
  } # for i 
  if (zz ==1) fileCompletenessID <- paste(text.add, ".completeExpression", sep="")
  comm.stack.pool <- paste(comm.stack.pool,  " && echo  >  ", fileCompletenessID, ".", zz, " & ", sep="")
} # for zz