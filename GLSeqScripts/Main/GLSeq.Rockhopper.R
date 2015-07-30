source("GLSeq.Util.R")
source("GLSeq.Alignment.Functions.R")
setwd(dest.dir)
# copy genome indices to the destimation dir: 
indCopy <- copyGenome(base.dir,rGenome,refFASTAname,dest.dir)
system(indCopy)
#
comm.stack.pool <- NULL # 




java -Xmx1200m -cp Rockhopper.jar Rockhopper -g <Reference Genome>  -o <OutputDirectory> -e false -SAM -TIME <FilePath> 

  put a % between paired files




for (zz in 1:nStreams) {
  # assembly and runing the system command, one library at a time:
  for (i in rangelist[[zz]]) {
    ###################
    # Alignment with SAM output
    ###################
    # Names of current fastq files:
    fq.left <- fqfiles.table[i,1]
    if (paired.end) fq.right <- fqfiles.table[i,2]
    name <- fqfiles.table[i,1]
    if (paired.end){
      name <- substr(name,1,length(name) - 5)
    } else{
      name <- substr(name,1,length(name) - 3)
    }
    this.resName <- paste(name, text.add, sep=".")
    unsorted.sam <- paste(this.resName, "unsorted", sep=".")
    #
  }
}