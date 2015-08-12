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
source("GLSeq.Alignment.Functions.R")
source("GLSeq.Util.R")

comm.stack.pool <- NULL
comm.stack.pools <- NULL
file.name.change <- NULL
####################################
# Copy genome indices to the destimation dir:
####################################
indCopy <- copy.genome(base.dir,rGenome,refFASTAname,dest.dir)
printOrExecute(indCopy,Condor)
####################

if (paired.end){
  for (zz in 1:nStreams) {
    for (i in rangelist[[zz]]){
      #
      # Tophat uses different file name conventions, so we need to switch them here.
      # TopHat wants paired files to be in the format
      # {FileNameBase}_1.fq while we normally put them in the format {FileNameBase}.1.fq
      # or
      # {FileNameBase}_2.fq while we normally put them in the format {FileNameBase}.2.fq
      #
      if (is.null(file.name.change)){
        file.name.change <- paste("mv",paste(dest.dir,fqfiles.table[i,1],sep=""))
      } else{
        file.name.change <- paste(file.name.change,"&& mv",paste(dest.dir,fqfiles.table[i,1],sep=""))
      }
      fqfiles.table[i,1] <- sub(".1.fq","_1.fq",fqfiles.table[i,1])
      file.name.change <- paste(file.name.change,paste(dest.dir,fqfiles.table[i,1],sep=""))
      file.name.change <- paste(file.name.change,"&& mv",paste(dest.dir,fqfiles.table[i,2],sep=""))
      fqfiles.table[i,2] <- sub(".2.fq","_2.fq",fqfiles.table[i,2])
      file.name.change <- paste(file.name.change,paste(dest.dir,fqfiles.table[i,2],sep=""))
    }
    # Change file names
    printOrExecute(file.name.change,Condor)
  }
}
# Index for TopHat Aligner with Bowtie2
#
# Prepare reference, what aligner to use, what reference to use, file name
#
#-f means FASTA input
####################
IndexOptions <- paste("-f")
index <- paste("bowtie2-build",IndexOptions,paste(dest.dir,refFASTAname,sep=""),paste(dest.dir,rGenome,sep=""))
printOrExecute(index,Condor)
#
for (zz in 1:nStreams) {
  comm.stack.pools <- NULL
  for (i in rangelist[[zz]]) {
    ###################
    # Alignment with SAM output
    ###################
    # Grabbing all the correct naming conventions and such.
    fq.left <- paste(dest.dir,fqfiles.table[i,1],sep="")
    if (paired.end) fq.right <- paste(dest.dir,fqfiles.table[i,2],sep="")
    name <- assign.name(fq.left,paired.end)
    this.resName <- assign.resName(name,text.add)
    #
    # -o Tells the aligner where to put the output, which is a unique folder based on the inputs.
    # This makes sure we don't overwrite any results and can pull them individually from their folder to
    # the main folder for further processing/display, but retain some of TopHat's excellent post-run notes
    #
    tophat.output.dir <- paste(this.resName,"TopHat.Alignment",sep=".")
    alignmentOptions <- paste("-p 8","-o",tophat.output.dir)
    align <- paste(TopHat.path,alignmentOptions,paste(dest.dir,rGenome,sep=""),fq.left)
    if (paired.end){
      align <- paste(align,fq.right)
    }
    # After the alignment, we'll want to move some files around and rename some files to line up more with our standard naming conventions.
    acceptedHitsName <- paste(this.resName,".accepted_hits.bam",sep="")
    tophat.output.dir <- trailDirCheck(tophat.output.dir)
    file.rename.and.move <- paste("mv",paste(tophat.output.dir,"accepted_hits.bam",sep=""),acceptedHitsName)
    bam.index <- paste("samtools index", acceptedHitsName)
    #
    ###################
    # Converting the BAM
    # to name-sorted  sam:
    ###################
    #
    countable.sam <- countable.sam.name(this.resName)
    paired.arg <- paste(this.resName, "paired", sep=".")
    paired.bam <- paste(this.resName, "paired.bam", sep=".")
    # Let's avoid pipe here to save RAM
    # Directly write the samtools output to a new countable file after sorting
    # -n in the samtools sort indicates that the document will be sorted by read name rather than chromosome
    # -h in the samtools view indicates a header will be printed
    if (Condor){
      countable.comm <- paste("samtools sort -@ 6 -m 32G -n", acceptedHitsName, paired.arg, "&&", "samtools view -h", paired.bam, ">", countable.sam)
    } else{
      countable.comm <- paste("samtools sort -n", acceptedHitsName, paired.arg, "&&", "samtools view -h", paired.bam, ">", countable.sam)
    }
    #
    count.comm <- ""
    if (counting == "counting"){
      setwd(base.dir)
      source("GLSeq.Counting.R")
    }
    #
    # Add any counting
    accepted.index <- paste(acceptedHitsName,".bai",sep="")
    remove.indexes <- paste(dest.dir,rGenome,".*.bt2",sep="")
    # Need to move the files that are made back into the normal folder
    if (is.null(comm.stack.pools)){
      comm.stack.pools <- paste(align)
    } else{
      comm.stack.pools <- paste(comm.stack.pools,'&&',align)
    }
    comm.stack.pools <- paste(comm.stack.pools,"&&",file.rename.and.move)
    comm.stack.pools <- paste(comm.stack.pools,"&&",bam.index)
    comm.stack.pools <- paste(comm.stack.pools,"&&",countable.comm)
    if (count.comm != "") comm.stack.pools <- paste(comm.stack.pools,"&&",count.comm)
    #
  }
  if (is.null(comm.stack.pool)){
    comm.stack.pool <- paste(comm.stack.pools,"&")
  } else{
    comm.stack.pool <- paste(comm.stack.pool,comm.stack.pools,"&")
  }
}
comm.stack.pool <- paste(comm.stack.pool,"wait")