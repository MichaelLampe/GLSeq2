source("GLSeq.Util.R")
source("GLSeq.Alignment.Functions.R")

hisat.directory <- trailDirCheck(hisat.directory)

comm.stack.pool <- NULL
comm.stack.pools <- NULL
indCopy <- copy.genome(base.dir,rGenome,refFASTAname,dest.dir)
printOrExecute(indCopy,Condor)

hisat.builder.location <- paste(hisat.directory,"hisat-build",sep="")
index <- paste("python",hisat.builder.location,paste(dest.dir,refFASTAname,sep=""),paste(dest.dir,rGenome,sep=""))
printOrExecute(index,Condor)

for (zz in 1:nStreams) {
  # assembly and runing the system command, one library at a time:
  for (i in rangelist[[zz]]) {
    ###################
    # Alignment with SAM output
    ###################
    # names of current fastq files:
    fq.left <- paste(dest.dir,fqfiles.table[i,1],sep="")
    if (paired.end) fq.right <- paste(dest.dir,fqfiles.table[i,2],sep="")
    name <- assign.name(fq.left,paired.end)
    this.resName <- assign.resName(name,text.add)
    countable.sam <- countable.sam.name(this.resName)
    ###################
    # Alignment
    ###################
    if (Condor){
      alignment.options <- paste("-p 12 -x",paste(dest.dir,rGenome,sep=""))
    } else{
      alignment.options <- paste("-x",paste(dest.dir,rGenome,sep=""))
    }
    hisat.location <- paste(hisat.directory,"hisat",sep="")
    if (paired.end){
      align <- paste("perl",hisat.location,alignment.options,"-1",fq.left,"-2",fq.right,"-S",countable.sam)
    } else{
      align <- paste("perl",hisat.location,alignment.options,"-U",fq.left,"-S",countable.sam)
    }
    ###################
    # Counting Step
    ###################
    count.comm <- ""
    if (counting == "counting"){
      setwd(base.dir)
      source("GLSeq.Counting.R")
    }
    ###################
    # Command Construction
    ###################
    comm.i <- paste(align)
    comm.i <- paste(comm.i, "&&", count.comm)
    # For the very first assembly in the stack:
    if (i == rangelist[[zz]][1])  comm.stack.pools <- paste(comm.i)
    # For subsequent assemblies of every stack:
    if (i != rangelist[[zz]][1])  comm.stack.pools <- paste(comm.stack.pools,"&&",comm.i)
    #
  }
  if (is.null(comm.stack.pool)){
    comm.stack.pool <- paste(comm.stack.pools,"&")
  } else{
    comm.stack.pool <- paste(comm.stack.pool,comm.stack.pools,"&")
  }
}
comm.stack.pool <- paste(comm.stack.pool,"wait")