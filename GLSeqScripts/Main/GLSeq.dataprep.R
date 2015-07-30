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
# Create QC Folder
qcFolder <-  create.QC.folder(dest.dir,text.add)

if (unzipped){
  fqFiles.unzip <- get.files.unzipped(raw.dir)
  fqFiles <- prepare.unzipped.file.names(fqFiles.unzip)
  copy.files.to.dest.unzipped(raw.dir,dest.dir)
}
if (!unzipped){
  fqFiles.zip <- get.files.zipped(raw.dir)
  fqFiles <- prepare.zipped.file.names(fqFiles.zip)
  copy.files.to.dest.zipped(raw.dir,dest.dir)
  for (i in 1:length(fqFiles.zip)){
    unzip.comm <- unzip.gz.files(fqFiles.zip[i])
    try(system(unzip.comm))
  }
}
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
  rangelist.Dataprep <- chunk.data.files.presplit(fqFiles,nStreamsDataPrep)
}
if (!presplit){
  rangelist.Dataprep <- chunk.data.files.unsplit(fqFiles,nStreamsDataPrep)
}

# Copy the artificial.fq file in
copy.artificial.fq(base.dir,artificial.fq,dest.dir)
# Construct command to be run.
for (zz in 1:nStreamsDataPrep) {
  if (zz==1) comm.pool <- "date"
  if (zz!=1) comm.pool <- paste(comm.pool,"date")
  for (j in rangelist.Dataprep[[zz]]) {
    #
    if (paired.end){
      #
      # Unsplit
      if (!presplit){
        # Add to command pool
        comm.pool <- paste(comm.pool,"&&",split.unsplit.files.PE(dest.dir,fqFiles[j]))
      }
      if (readTrim){
        # The trim command
        trimCommand <- trimAssemble.PE(fqFiles[j], trimPath, qScores, trimhead, artificial.fq,trimMin)
        comm.pool <- paste(comm.pool,"&&",trimCommand)
        # Modifies the names of a few files to retain naming
        file.shuffle <- file.shuffle.PE(fqFiles[j])
        comm.pool <- paste(comm.pool,"&&",file.shuffle)
        # Quality control check via Fastqc of dirty file
        preQC <- preQualityCheck.PE(fastqcPath,fqFiles[j],qcFolder)
        comm.pool <- paste(comm.pool,"&&",preQC)
        remove.command <- remove.unneeded.files(fqFiles[j])
        comm.pool <- paste(comm.pool,"&&",remove.command)
        move.file.log <- move.paired.files.PE((fqFiles[j]),qcFolder)
        comm.pool <- paste(comm.pool,"&&",move.file.log)
      }   
      # Quality control check via Fastqc of result files
      postQC <- postQualityCheck.PE(fastqcPath,fqFiles[j],qcFolder)
      comm.pool <- paste(comm.pool,"&&",postQC)
    }
    #
    if (!paired.end){
      if (readTrim){
        # The trim command
        trimCommand <- trimAssemble.SE(fqFiles[j], trimPath, qScores, trimhead, artificial.fq,trimMin)
        comm.pool <- paste(comm.pool,"&&",trimCommand)
        # Modifies the names of a few files to retain naming
        file.shuffle <- file.shuffle.SE(fqFiles[j])
        comm.pool <- paste(comm.pool,"&&",file.shuffle)
        # Quality control check via Fastqc of dirty file
        preQC <- preQualityCheck.SE(fastqcPath,fqFiles[j],qcFolder)
        comm.pool <- paste(comm.pool,"&&",preQC)
        remove.command <- remove.unneeded.files.SE(fqFiles[j])
        comm.pool <- paste(comm.pool,"&&",remove.command)
        move.file.log <- move.paired.files.SE((fqFiles[j]),qcFolder)
        comm.pool <- paste(comm.pool,"&&",move.file.log)
                                        
      }
      # Quality control check via Fastqc of result files
      postQC <- postQualityCheck.SE(fastqcPath,fqFiles[j],qcFolder)
      comm.pool <- paste(comm.pool,"&&",postQC)
    }
  }
  comm.pool <- paste(comm.pool,"&")
}
store.artificial <- store.artificial.seqs.file(artificial.fq,qcFolder)
comm.pool <- paste(comm.pool,"wait","&&", store.artificial)
try(system(comm.pool))