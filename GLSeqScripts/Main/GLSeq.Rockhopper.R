source("GLSeq.Util.R")
source("GLSeq.Alignment.Functions.R")
# copy genome indices to the destimation dir:
if (seqPlatform != "illumina"){
  stop("Only illumina platform is supported for Rockhopper currently.")
}
#

rGenome <- trailDirCheck(rGenome)
ref.dir <- paste(base.dir, rGenome, sep="")
rGenomeDestination <- paste(dest.dir,"ReferenceGenome",sep="")
indCopy <- paste("mkdir",rGenomeDestination,"&&","cp",paste(ref.dir,refFASTAname,sep=""),rGenomeDestination)
printOrExecute(indCopy,Condor)
#
# Small script I wrote that breaks fasta files into individual FNA files in the format that Rockhopper likes.
convertReferenceGenome <- paste("python",paste(base.dir,"rhFnaConverter.py",sep=""),rGenomeDestination,refFASTAname,dest.dir)
printOrExecute(convertReferenceGenome,Condor)
# Creates ptt files from a single gtf file
ref.dir <- paste(base.dir, rGenome, sep="")
ref.dir <- trailDirCheck(ref.dir)
refGtf <- paste(ref.dir,refGFFname,sep="")
refCopy <- paste("cp ",refGtf," ",dest.dir, sep="")
printOrExecute(refCopy,Condor)
gtfToPtt <- paste("python",paste(base.dir,"GtfToPtt.py",sep=""),refGtf,rGenomeDestination)
printOrExecute(gtfToPtt,Condor)
# check the environment over.  If a GTF file has 0 transcripts for a certain FASTA sequence, rockhopper will crash if that folder remains.
checkEnvironment <- paste("python checkRhEnvironment.py",rGenomeDestination)
printOrExecute(checkEnvironment,Condor)
# This needs to catch the output doesnt it
#
comm.stack.pool <- NULL #
comm.stack.pools <- NULL
# Checks the rockhopper path
Rockhopper.path <- trailDirCheck(Rockhopper.path)
for (zz in 1:nStreams) {
  comm.stack.pool <- NULL
  # assembly and runing the system command, one library at a time:
  for (i in rangelist[[zz]]) {
    print("FQ FILES")
    print(fqfiles.table[i,1])
    ###################
    # Alignment with SAM output
    ###################
    # Names of current fq files:
    # Need to be renamed to fastq if not already
    if (grepl(".fq$",fqfiles.table[i,1])){
      fastq.name <- substr(fqfiles.table[i,1],1,nchar(fqfiles.table[i,1]) - 2)
      fastq.name <- paste(fastq.name,"fastq",sep="")
      rename.fq.file <- paste("mv",paste(dest.dir,fqfiles.table[i,1],sep=""),paste(dest.dir,fastq.name,sep=""))
      fqfiles.table[i,1] <- fastq.name
      printOrExecute(rename.fq.file,Condor)
    }
    # Paired ended need both adjusted
    if(paired.end){
      if (grepl(".fq$",fqfiles.table[i,2])){
        fastq.name <- substr(fqfiles.table[i,2],1,nchar(fqfiles.table[i,2]) - 2)
        fastq.name <- paste(fastq.name,"fastq",sep="")
        rename.fq.file <- paste("mv",paste(dest.dir,fqfiles.table[i,2],sep=""),paste(dest.dir,fastq.name,sep=""))
        fqfiles.table[i,2] <- fastq.name
        printOrExecute(rename.fq.file,Condor)
      }
    }
    fq.left <- fqfiles.table[i,1]
    if (paired.end) fq.right <- fqfiles.table[i,2]
    name <- fqfiles.table[i,1]
    if (paired.end){
      name <- substr(name,1,nchar(name) - 8)
    } else{
      name <- substr(name,1,nchar(name) - 6)
    }
    this.resName <- paste(paste(dest.dir,name,sep=""), text.add, sep=".")
    # The reference genomes are the folders within rGenomeDestination.
    # We'll need to use a way to at runtime get the folders within that destination after the rhFnaConverter script has run.
    # The solution I came up with was to just run a python script that could determine all of this stuff at runtime and then call
    if (paired.end){
      # Rockhopper denotes paired files by the % between them
      rockhopper.files <- paste(dest.dir,fq.left,"%",dest.dir,fq.right,sep="")
    } else{
      # Nonpaired files are just run as normal.
      rockhopper.files <- paste(dest.dir,fq.left,sep="")
    }
    rockhopper.output.folder <- paste(this.resName,"RockhopperResults/",sep=".")

    if (paired.end){
      if (is.null(libstrand)) {
        strand.command <- "-s false"
      }
      else if (libstrand == "F") {
        strand.command <- "-fr -s true"
      }
      else if (libstrand == "R") {
        strand.command <- "-rf -s true"
      }
    } else{
      if (is.null(libstrand)) {
        strand.command <- "-s false"
      }
      else if (libstrand == "F") {
        strand.command <- "-c false"
      }
      else if (libstrand == "R") {
        strand.command <- "-c true"
      }
    }

    rockhopper.align <- paste("python",paste(base.dir,"rockhopperWrapper.py",sep=""),Rockhopper.path,dest.dir,rockhopper.files,rockhopper.output.folder,strand.command)
    countable.sam <- countable.sam.name(this.resName)
    create.countable <- paste("mv",paste(rockhopper.output.folder,"dirty.sam",sep=""),countable.sam)

    count.comm <- ""
    if (counting == "counting"){
      setwd(base.dir)
      source("GLSeq.Counting.R")
    }

    # If the comm stack is null,just add rockhopper, otherwise add to the list via &&
    if (is.null(comm.stack.pool)){
      comm.stack.pool <- paste(rockhopper.align,"&&",create.countable)
    } else{
      comm.stack.pool <- paste(comm.stack.pool,"&&",rockhopper.align,"&&",create.countable)
    }
    if (count.comm != ""){
      comm.stack.pool <- paste(comm.stack.pool,"&&",count.comm)
    }
  }
  # Add it to a unique background pool
  if (is.null(comm.stack.pools)){
    comm.stack.pools <- paste(comm.stack.pool,"&")
  } else{
    comm.stack.pools <- paste(comm.stack.pools,comm.stack.pool,"&")
  }
}

# This just retains the consistency of the program
comm.stack.pool <- paste(comm.stack.pools,"wait")