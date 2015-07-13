#########################################################
# Great Lakes Seq package for low-level processing of RNA-Seq data
# Trimmomatic-BWA-HTSeq quantification of the expression values
# TopHat Alignment Method
#
# Date: June 2015
# Author: Michael Lampe
#########################################################
#
#
source("GLSeq.Util.R")
setwd(dest.dir)
#
file.name.change <- "date"
if (paired.end){
  for (zz in 1:nStreams) {
    for (i in rangelist[[zz]]){
      #
      # Tophat uses different file name conventions, so we need to switch them here.
      #
      file.name.change <- paste(file.name.change,"&&","mv",fqfiles.table[i,1])
      fqfiles.table[i,1] <- sub(".1.fq","_1.fq",fqfiles.table[i,1])
      file.name.change <- paste(file.name.change,fqfiles.table[i,1])
      file.name.change <- paste(file.name.change,"&&","mv",fqfiles.table[i,2])
      fqfiles.table[i,2] <- sub(".2.fq","_2.fq",fqfiles.table[i,2])
      file.name.change <- paste(file.name.change,fqfiles.table[i,2])
    }
    # Change file names
    system(file.name.change)
  }
}
#
####################################
# Copy genome indices to the destimation dir: 
####################################
ref.dir <- paste(base.dir, rGenome, sep="")
indCopy <- paste("cd ", ref.dir, " && cp ",refFASTAname," ",dest.dir, sep="")
system(indCopy)
####################
# Index for TopHat Aligner with Bowtie2
#
# Prepare reference, what aligner to use, what reference to use, file name
#
#-f means FASTA input
####################
IndexOptions <- paste("-f")
index <- paste("bowtie2-build",IndexOptions,refFASTAname,rGenome)
system(index)
#
#
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
    ##
    # -o Tells the aligner where to put the output, which is a unique folder based on the inputs.  
    # This makes sure we don't overwrite any results and can pull them individually from their folder to
    # the main folder for further processing/display, but retain some of TopHat's excellent post-run notes
    #
    tophat.output.dir <- paste(this.resName,"TopHat.Alignment",sep=".")
    alignmentOptions <- paste("-o",tophat.output.dir)
    align <- paste(TopHat.path,alignmentOptions,rGenome,fq.left)
    if (paired.end){
      align <- paste(align,fq.right)
    }
    acceptedHitsName <- paste("Accepted_Hits.",this.resName,".bam",sep="")
    file.rename.and.move <- paste("cd",tophat.output.dir,"&&","mv","accepted_hits.bam",acceptedHitsName,"&&","mv",acceptedHitsName,dest.dir,"&&","cd",dest.dir)
    
    bam.index <- paste("samtools index", acceptedHitsName)
    #
    ###################
    # Converting the BAM
    # to name-sorted  sam:  
    ###################
    #
    countable.sam <- paste(this.resName, "countable.sam", sep=".")
    paired.arg <- paste(this.resName, "paired", sep=".")
    paired.bam <- paste(this.resName, "paired.bam", sep=".")
    # let's avoid pipe here to save RAM
    # Directly write the samtools output to a new countable file after sorting
    # -n in the samtools sort indicates that the document will be sorted by read name rather than chromosome
    # -h in the samtools view indicates a header will be printed
    countable.comm <- paste("samtools sort -n", acceptedHitsName, paired.arg, "&&", "samtools view -h", paired.bam, ">", countable.sam)
    #
    count.comm <- ""
    if (counting == "counting"){
      setwd(base.dir)
      source("GLSeq.Counting.R")
    }
    #
    # Add any counting
    # Need to move the files that are made back into the normal folder
    if (i == 1) comm.stack.pool <- paste("cd",dest.dir) 
    if (i != 1) comm.stack.pool <- paste(comm.stack.pool,"&&","cd",dest.dir)
    comm.stack.pool <- paste(comm.stack.pool,'&&',align)
    comm.stack.pool <- paste(comm.stack.pool,"&&",file.rename.and.move)
    comm.stack.pool <- paste(comm.stack.pool,"&&",bam.index)
    comm.stack.pool <- paste(comm.stack.pool,"&&",countable.comm)
    if (count.comm != "") comm.stack.pool <- paste(comm.stack.pool,"&&",count.comm)
    #    
  }
  comm.stack.pool <- paste(comm.stack.pool,"&")
}
comm.stack.pool <- paste(comm.stack.pool,"wait")