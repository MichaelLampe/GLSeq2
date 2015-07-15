#########################################################
# Great Lakes Seq package for low-level processing of RNA-Seq data
# Oleg Moskvin; info@scienceforever.com 
# April 2013 
#########################################################
# 
# Preparing the raw data
#
#########################################################
#
###########################
# ID of the computation run (text.add object) should be 
# supplied to this script directly (by GLSeq.top.R) as 
# the first and only argument
###########################
#
args <- commandArgs(trailingOnly = TRUE)
text.add <- as.character(args[1])
dest.dir <- as.character(args[2])
attrPath <- as.character(args[3])
source(attrPath)
source("GLSeq.Dataprep.Functions.R")
#
###########################
setwd(dest.dir)
###########################
# Files indicating readiness of the libraries (i.e. split FQ files):
###########################
files2watch.dataprep <- NULL
fqFiles.zip <- NULL
fqFiles.unzip <- NULL
###########################
load.dataFile(text.add)

# Grab the files
if (unzipped){
  fqFiles.unzip <- get.files.unzipped(raw.dir)
  fqFiles <- prepare.unzipped.file.names(fqFiles.unzip)
}
if (!unzipped){
  fqFiles.zip <- get.files.zipped(raw.dir)
  fqFiles <- prepare.zipped.file.names(fqFiles.zip)
}
#
#
#
check.ifFiles(fqFiles.zip,fqFiles.unzip)
#
# Check the values of two variables.
presplit <- check.presplit(paired.end,presplit)
nStreamsDataPrep <- check.nStreamsDataPrep(fqFiles,nStreamsDataPrep)
#
# Chunk the files into segments to be run in parallel.
#
# Note: This is a safe operation, nStreamsDataPrep can exceed the number
# of files and not cause a crash (It is checked above and within the function so it is modular, but explicit also)
if (presplit){
  rangelist.Dataprep <- chunk.data.files.presplit(fqFiles,nStreamDataPrep)
}
if (!presplit){
  rangelist.Dataprep <- chunk.data.files.unsplit(fqFiles.nStreamsDataPrep)
}

# Copy the artificial.fq file in
copy.artificial.fq(base.dir,artificial.fq,dest.dir)
# Construct command to be run.
for (zz in 1:nStreamsDataPrep) {
  if (zz==1) comm.pool <- "date"
  if (zz!=1) comm.pool <- paste(comm.pool,"date")
  for (j in rangelist.Dataprep[[zz]]) {
    #
    # GENERAL PROCESS
    #
    # Constructs a command to copy over a file.
    if (unzipped){
      copy.comm <- paste(copy.comm,"&&",copy.file(copy.comm,raw.dir,fqFiles.unzip[j],dest.dir))
      #
      if(presplit){
        # Makes sure both are copied over in the same command so there is not a parallelization issue here.
        copy.comm <- paste(copy.comm,"&&",copy.file(copy.comm,raw.dir,fqFiles.unzip[j+1],dest.dir))
      }
    }
    if (!unzipped){
      # Constructs a command to unzip .gz zipped files.
      copy.comm <- paste(copy.comm,"&&",copy.file(copy.comm,raw.dir,fqFiles.zip[j],dest.dir))
      copy.comm <- paste(copy.comm,"&&",unzip.gz.files(fqFile.zip[j]))
      #
      if(presplit){
        # Makes sure both are copied over in the same command so there is not a parallelization issue here.
        copy.comm <- paste(copy.comm,"&&",copy.file(copy.comm,raw.dir,fqFiles.zip[j+1],dest.dir))
        copy.comm <- paste(copy.comm,"&&",unzip.gz.files(fqFile.zip[j+1]))
      }
    }
    comm.pool <- paste(comm.pool,"&&",copy.comm)
    #
    # Paired Ended Files
    #
    if (paired.end){
      #
      # Unsplit
      if (!presplit){
        # Add to command pool
        comm.pool <- paste(comm.pool,"&&",split.unsplit.files.PE(dest.dir,fqFile[j]))
      }
      if (readTrim){
        # The trim command
        trimCommand <- trimAssemble.PE(fqFile[j], trimPath, qScores, trimhead, artificial.fq)
        comm.pool <- paste(comm.pool,"&&",trimCommand)
        # Modifies the names of a few files to retain naming
        file.shuffle <- file.shuffle.PE(fqFile[j])
        comm.pool <- paste(comm.pool,"&&",file.shuffle)
        # Quality control check via Fastqc of dirty file
        preQC <- preQualityCheck.PE(fastqcPath,fqFile[j])
        comm.pool <- paste(comm.pool,"&&",preQC)
        # Quality control check via Fastqc of result files
        postQC <- postQualityCheck.PE(fastqcPath,fqFile[j])
        comm.pool <- paste(comm.pool,"&&",postQC) 
      }
      if (!readTrim){
        postQC <- postQualityCheck.PE(fastqcPath,fqFile[j])
        comm.pool <- paste(comm.pool,"&&",postQC) 
      }
      if(presplit){
        # Move onto the next set of pairs if presplit data.
        j = j + 1
      }
    }
    if (!paired.end){
      if (readTrim){
        # The trim command
        trimCommand <- trimAssemble.SE(fqFile[j], trimPath, qScores, trimhead, artificial.fq)
        comm.pool <- paste(comm.pool,"&&",trimCommand)
        # Modifies the names of a few files to retain naming
        file.shuffle <- file.shuffle.PE(fqFile[j])
        comm.pool <- paste(comm.pool,"&&",file.shuffle)
        # Quality control check via Fastqc of dirty file
        preQC <- preQualityCheck.PE(fastqcPath,fqFile[j])
        comm.pool <- paste(comm.pool,"&&",preQC)
        # Quality control check via Fastqc of result files
        postQC <- postQualityCheck.PE(fastqcPath,fqFile[j])
        comm.pool <- paste(comm.pool,"&&",postQC) 
      }
      if (!readTrim){
        postQC <- postQualityCheck.PE(fastqcPath,fqFile[j])
        comm.pool <- paste(comm.pool,"&&",postQC) 
      }
    }
  }
  comm.pool <- paste(comm.pool,"&")
}
dataReady.signal <- paste("echo > ", text.add, ".DataReady", sep="") 
comm.pool <- paste(comm.pool,"wait","&&",dataReady.signal)
system(comm.pool)