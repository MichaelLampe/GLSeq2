source("GLSeq.Util.R")
source("GLSeq.Alignment.Functions.R")

# The final return pool for the program
comm.stack.pool <- NULL
comm.stack.pools <- NULL
indCopy <- copy.genome(base.dir,rGenome,refFASTAname,dest.dir)
printOrExecute(indCopy,Condor)
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
  index <- paste("bowtie-build",IndexOptions,paste(dest.dir,refFASTAname,sep=""),paste(dest.dir,rGenome,sep=""))
}
if (aAlgor == "Bowtie2"){
  index <- paste("bowtie2-build",IndexOptions,paste(dest.dir,refFASTAname,sep=""),paste(dest.dir,rGenome,sep=""))
}
printOrExecute(index,Condor)
#
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
    # --sam gives a sam file output, -t gives time information for the process.
    # --chunkmbs may fix paired end read problems skipping reads due to memory exhaustion
    # Default for chunkmbs = 64, 512 is 8x larger.  No known problems created by this option.
    ###################
    if (aAlgor == "Bowtie"){
      if (Condor){
        BowtieOptions <- paste("-S","-p 4","-t","-q","--chunkmbs","512",alignmentSpecialOptions)
      } else{
        BowtieOptions <- paste("-S","-t","-q","--chunkmbs","512",alignmentSpecialOptions)
      }
      if (paired.end){
        align <- paste("bowtie",BowtieOptions,paste(dest.dir,rGenome,sep=""),"-1",fq.left,"-2",fq.right,countable.sam)
      }else{
        align <- paste("bowtie",BowtieOptions,paste(dest.dir,rGenome,sep=""),fq.left,countable.sam)
      }
    }
    ###################
    if (aAlgor == "Bowtie2"){
      # All these options are necessary just to get it to conform with the RSEM counting protocol... (Except -t and -q)
      # Found on the RSEM website (http://deweylab.biostat.wisc.edu/rsem/rsem-calculate-expression.html)
      # This is under the --bowtie2 option that is coupled with the RSEM package.
      # Might be worth moving off the RSEM counting protocol if we believe Bowtie will better serve us with the current
      # counting protocols.
      if (Condor){
        # I've set Bowtie to have 4 CPUS when using condor, so let's exploit that.
        # This is under the assumption that any leftover cores would be able to be used by other processes while running in parallel.
        # Thus, we'll be able to fill up the servers as efficiently as posible with the LCD, which is 4.
        # Based on https://www.biostars.org/p/92366/ as well showing a non-linear speedup of runs. (4 -> 8 is not half time)
        Bowtie2Options <- paste("-t","-p 4","-q", alignmentSpecialOptions)
      } else{
        Bowtie2Options <- paste("-t","-q", alignmentSpecialOptions)
      }
      if (paired.end){
        # These are options normally used by the RSEM peeps, see the same link as above for documentation on this.
        Bowtie2Options <- paste(Bowtie2Options,"--no-mixed","--no-discordant")
        align <- paste("bowtie2",Bowtie2Options,paste(dest.dir,rGenome,sep=""),"-1",fq.left,"-2",fq.right,"-S",countable.sam)
      } else{
        align <- paste("bowtie2",Bowtie2Options,paste(dest.dir,rGenome,sep=""),"-U",fq.left,"-S",countable.sam)
      }
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