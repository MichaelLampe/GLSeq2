source("GLSeq.Util.R")
source("GLSeq.Alignment.Functions.R")

# The final return pool for the program
comm.stack.pool <- NULL
comm.stack.pools <- NULL
# Makes the reference genome directory
rGenome <- trailDirCheck(rGenome)
reference.genome.folder <- paste(dest.dir,rGenome,sep="")
create.directory <- paste("mkdir",reference.genome.folder)
printOrExecute(create.directory,Condor)
# Copies the fasta file into reference genome
reference.fasta.location <- paste(base.dir,rGenome,refFASTAname,sep="")
indCopy <- paste("cp",reference.fasta.location,reference.genome.folder)
printOrExecute(indCopy,Condor)

# Index
star.path <- "/home/GLBRCORG/mrlampe/STAR/STAR"

index <- paste(star.path,"--runMode genomeGenerate --genomeDir",reference.genome.folder,"--genomeFastaFiles",paste(reference.genome.folder,refFASTAname,sep=""))
if (Condor){
  paste(index,"--runThreadN 6")
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
    star.output.dir <- paste(this.resName,"STAR.Alignment/",sep=".")
    printOrExecute(paste("mkdir",star.output.dir),Condor)
    if (paired.end){
      alignment.command <- paste(star.path,"--genomeDir",reference.genome.folder,alignmentSpecialOptions,"--readFilesIn",paste(fq.left,",",fq.right,sep=""),"--outFileNamePrefix",star.output.dir)
    } else{
      alignment.command <- paste(star.path,"--genomeDir",reference.genome.folder,alignmentSpecialOptions,"--readFilesIn",fq.left,"--outFileNamePrefix",star.output.dir)
    }
    # 6 Cores, roughly 24GB of ram they estimate for this.
    if (Condor) alignment.command <- paste(alignment.command,"--runThreadN 6")
    generate.countable.sam <- paste("mv",paste(star.output.dir,"Aligned.out.sam",sep=""),countable.sam)
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
    comm.i <- paste(alignment.command, "&&", generate.countable.sam)
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