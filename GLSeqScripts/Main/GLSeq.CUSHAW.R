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
source("GLSeq.Alignment.Functions.R")

comm.stack.pool <- NULL

indCopy <- copy.genome(base.dir,rGenome,refFASTAname,dest.dir)
printOrExecute(indCopy,Condor)
# This will index the CUSHAW with your FASTA file, necessary.
# Previously the user needed to do this, but it should be easier if we just take care of
# This stuff on our end during each run (It is fast && cheap anyway)
index <- paste(CUSHAW.index.path,"-a","bwtsw",paste(dest.dir,refFASTAname,sep=""))
####################################
# ~Hot~ system command to index aligner
####################################
printOrExecute(index,Condor)
####################################
### Let's speed this up.
### If using the GPU, alignment occurs in sequence so as to not overload the GPU memory (Each process is ~6GB)
### Our NVidia Titan can sometimes run two at a time without alignment corruption, but this is a safer route.
####################################
sam.create <- NULL
if (aAlgor == "Cushaw_GPU"){
  for (zz in 1:nStreams) {
    for (i in rangelist[[zz]]) {
      ###################
      # Alignment with SAM output
      ###################
      # names of current fastq files:
      fq.left <- paste(dest.dir,fqfiles.table[i,1],sep="")
      if (paired.end) fq.right <- paste(dest.dir,fqfiles.table[i,2],sep="")
      name <- paste(dest.dir,fqfiles.table[i,1],sep="")
      if (paired.end){
        name <- substr(name,1,nchar(name) - 5)
      } else{
        name <- substr(name,1,nchar(name) - 3)
      }
      this.resName <- paste(name, text.add, sep=".")
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
      if (paired.end)  create <- paste(CUSHAW.GPU.path, "-r", paste(dest.dir,refFASTAname,sep=""), "-q", fq.right, fq.left, "-o", unsorted.sam,"-t", nCores)
      if (!(paired.end)) create <- paste(CUSHAW.GPU.path, "-r", paste(dest.dir,refFASTAname,sep=""), "-f", fq.left, "-o", unsorted.sam,"-t", nCores)
      # Checks to make sure that this process actually is GPU
      # The is.null is there for later on implementation of parallel GPU runs. (Hopefully!)
      if (is.null(sam.create)){
        sam.create <- paste(create)
      } else{
        sam.create <- paste(sam.create,"&&",create)
      }
      # Fixes the order of processing so that
      # The SAM file gets created when the user
      # Preps data, but it still executes in sequence
      #
    } # i
  } #nStreams (zz)
  printOrExecute(sam.create,Condor)
} #GPU.Accel
comm.stack.pool <- NULL
comm.stack.pools <- NULL
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
    fq.left <- paste(dest.dir,fqfiles.table[i,1],sep="")
    if (paired.end) fq.right <- paste(dest.dir,fqfiles.table[i,2],sep="")
    name <- paste(dest.dir,fqfiles.table[i,1],sep="")
    if (paired.end){
      name <- substr(name,1,nchar(name) - 5)
    } else{
      name <- substr(name,1,nchar(name) - 3)
    }
    this.resName <- paste(name, text.add, sep=".")
    unsorted.sam <- paste(this.resName, "unsorted", sep=".")
    #
    ###################
    # If no GPU, we can run all of the above files in parallel
    ###################
    # Condor is set to use 8 cores
    if (Condor){
      nCores <- 8
    }
    if (aAlgor == "Cushaw"){
      if (paired.end) sam.create<- paste(CUSHAW.path, "-r", paste(dest.dir,refFASTAname,sep=""), "-q", fq.left, fq.right, "-o", unsorted.sam, "-t", nCores)
      if (!paired.end) sam.create <- paste(CUSHAW.path, "-r", paste(dest.dir,refFASTAname,sep=""), "-f", fq.left, "-o", unsorted.sam,"-t", nCores)
    }
    ###################
    # Convert BAM to sorted BAM file
    ###################
    sorted <- paste(this.resName,"sorted",sep=".")
    unsorted.bam <- paste(this.resName,"unsorted.bam",sep=".")
    sorted <- paste(this.resName, "sorted", sep=".")
    bam.create <- paste("samtools view",unsorted.sam,"-uS -o",unsorted.bam)
    if (Condor){
      bam.sort <- paste("samtools sort -@ 6 -m 32G -n ",unsorted.bam,sorted)
    } else{
      bam.sort <- paste("samtools sort -n ",unsorted.bam,sorted)
    }
    bam.index <- paste("samtools index", sorted) # System command #7
    ###################
    # Convert back to SAM file
    ###################
    countable.sam <- paste(this.resName,"countable","sam",sep=".")
    # This removes certain read errors that can cause counting programs to crash
    # Is the argument to an AWK command that the data is piped into before being written.
    # This is only needed when using paired.ended data, so you can see we avoid that pipe when
    # we are utilizing SE data.
    # Shell script made from a comment found at https://www.biostars.org/p/108702/
    sorted.bam <- paste(this.resName,"sorted.bam",sep=".")
    unclean.sam <- paste(this.resName,"unclean.sam",sep=".")
    if(paired.end) {
      convert.to.sam <- paste ("samtools view -h",sorted.bam,">",unclean.sam)
      filter.sam <- paste("bash",paste(base.dir,"CushawCorrection.sh",sep=""),unclean.sam,countable.sam)
    }
    if (!paired.end) convert.to.sam <- paste("samtools view",sorted.bam,"-h -o ",countable.sam)
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
    if (aAlgor == "Cushaw") comm.i <- paste(sam.create,"&&",bam.create,"&&",bam.sort,"&&",bam.index,"&&",convert.to.sam)
    if (aAlgor == "Cushaw_GPU") comm.i <- paste(bam.create,"&&",bam.sort,"&&",bam.index,"&&",convert.to.sam)
    if (paired.end) comm.i <- paste(comm.i,"&&",filter.sam)
    #
    if (count.comm != "") comm.i <- paste(comm.i, "&&", count.comm)
    #
    # For the very first assembly in the stack (i = 1)
    if (i == rangelist[[zz]][1])  comm.stack.pool <- paste(comm.i)
    #
    # For subsequent assemblies of every stack (i > 1)
    if (i != rangelist[[zz]][1])  comm.stack.pool <- paste(comm.stack.pool,"&&", comm.i)
  } # for i

  if (is.null(comm.stack.pools)){
    comm.stack.pools <- paste(comm.stack.pool,"&")
  } else{
    comm.stack.pools <- paste(comm.stack.pools,comm.stack.pool,"&")
  }
}
comm.stack.pool <- paste(comm.stack.pools,"wait")