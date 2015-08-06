source("GLSeq.Util.R")
source("GLSeq.Alignment.Functions.R")
setwd(dest.dir)
# copy genome indices to the destimation dir:

rGenome <- trailDirCheck(rGenome)
ref.dir <- paste(base.dir, rGenome, sep="")
rGenomeDestination <- paste(dest.dir,"ReferenceGenome",sep="")
indCopy <- paste("mkdir",rGenomeDestination,"&&","cp",paste(ref.dir,refFASTAname,sep=""),rGenomeDestination)
printOrExecute(indCopy,Condor)
#
# Small script I wrote that breaks fasta files into individual FNA files in the format that Rockhopper likes.
convertReferenceGenome <- paste("python",paste(base.dir,"rhFnaConverter.py",sep=""),rGenomeDestination,refFASTAname)
printOrExecute(convertReferenceGenome,Condor)
# This needs to catch the output doesnt it
#

#
comm.stack.pool <- NULL #
comm.stack.pools <- NULL
# Checks the rockhopper path
Rockhopper.path <- trailDirCheck(Rockhopper.path)
for (zz in 1:nStreams) {
  comm.stack.pool <- NULL
  # assembly and runing the system command, one library at a time:
  for (i in rangelist[[zz]]) {
    ###################
    # Alignment with SAM output
    ###################
    # Names of current fq files:
    # Need to be renamed to fastq if not already
    print("HI")
    if (grepl(".fq$",fqfiles.table[i,1])){
      fastq.name <- substr(fqfiles.table[i,1],1,nchar(fqfiles.table[i,1]) - 2)
      fastq.name <- paste(fastq.name,"fastq",sep="")
      rename.fq.file <- paste("mv",fqfiles.table[i,1],fastq.name)
      fqfiles.table[i,1] <- fastq.name
      printOrExecute(rename.fq.file,Condor)
    }
    # Paired ended need both adjusted
    if(paired.end){
      if (grepl(".fq$",fqfiles.table[i,2])){
        fastq.name <- substr(fqfiles.table[i,2],1,nchar(fqfiles.table[i,2]) - 2)
        fastq.name <- paste(fastq.name,"fastq",sep="")
        rename.fq.file <- paste("mv",fqfiles.table[i,2],fastq.name)
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
    rockhopper.align <- paste("python",paste(base.dir,"rockhopperWrapper.py",sep=""),Rockhopper.path,dest.dir,rockhopper.files,paste(this.resName,"RockhopperResults",sep="."))
    # We currently have no use for the SAM file generated (I can't figure out how to count it..., so I'll just leave that stuff aside.)

    # If the comm stack is null,just add rockhopper, otherwise add to the list via &&
    if (is.null(comm.stack.pool)){
      comm.stack.pool <- rockhopper.align
    } else{
      comm.stack.pool <- paste(comm.stack.pool,"&&",rockhopper.align)
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
comm.stack.pool <- comm.stack.pools