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
source("GLSeq.Util.R")
source(attrPath)
source("GLSeq.Dataprep.Functions.R")
#
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
  unzip.comm <- NULL
  for (i in 1:length(fqFiles.zip)){
    if (is.null(unzip.comm)){
    unzip.comm <- unzip.gz.files(paste(dest.dir,fqFiles.zip[i],sep=""))
    } else{
      unzip.comm <- paste(unzip.comm,unzip.gz.files(paste(dest.dir,fqFiles.zip[i],sep="")))
    }
  }
  # Parallel Gunzip
  printOrExecute(unzip.comm,Condor)
}
# Most of the alignment stuffs adds the destination and such on their own, so let's not mess that up and return
# Just the relative file name
relative.fqFiles <- fqFiles
fqFiles <- convert.to.absolute.paths(dest.dir,fqFiles)
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
if (readTrim){
  copy.artificial.fq(base.dir,artificial.fq,dest.dir,Condor)
}
# Construct command to be run.
comm.pool <- NULL
for (zz in 1:nStreamsDataPrep) {
  comm.pools <- NULL
  for (j in rangelist.Dataprep[[zz]]) {
    #
    if (paired.end){
      #
      # Unsplit
      if (!presplit){
        # Add to command pool
        if (is.null(comm.pools)){
          comm.pools <- paste(split.unsplit.files.PE(base.dir,fqFiles[j]))
        }else{
          comm.pools <- paste(comm.pools,"&&",split.unsplit.files.PE(dest.dir,fqFiles[j]))
        }
        }
      if (readTrim){
        # The trim command
        trimCommand <- trimAssemble.PE(fqFiles[j], trimPath, qScores, trimhead, artificial.fq,trimMin)
        if (is.null(comm.pools)){
          comm.pools <- paste(trimCommand)
        } else{
          comm.pools <- paste(comm.pools,"&&",trimCommand)
        }
        # Modifies the names of a few files to retain naming
        file.shuffle <- file.shuffle.PE(fqFiles[j])
        comm.pools <- paste(comm.pools,"&&",file.shuffle)
        # Quality control check via Fastqc of dirty file
        preQC <- preQualityCheck.PE(fastqcPath,fqFiles[j],qcFolder)
        comm.pools <- paste(comm.pools,"&&",preQC)
        move.file.log <- move.paired.files.PE((fqFiles[j]),qcFolder)
        comm.pools <- paste(comm.pools,"&&",move.file.log)
      }
      # Quality control check via Fastqc of result files
      postQC <- postQualityCheck.PE(fastqcPath,fqFiles[j],qcFolder)
      if (is.null(comm.pools)){
        comm.pools <- paste(postQC)
      } else{
        comm.pools <- paste(comm.pools,"&&",postQC)
      }
    }
    #
    if (!paired.end){
      if (readTrim){
        # The trim command
        trimCommand <- trimAssemble.SE(fqFiles[j], trimPath, qScores, trimhead, artificial.fq,trimMin)
        if (is.null(comm.pools)){
          comm.pools <- paste(trimCommand)
        } else{
          comm.pools <- paste(comm.pools,"&&",trimCommand)
        }
        # Modifies the names of a few files to retain naming
        file.shuffle <- file.shuffle.SE(fqFiles[j])
        comm.pools <- paste(comm.pools,"&&",file.shuffle)
        # Quality control check via Fastqc of dirty file
        preQC <- preQualityCheck.SE(fastqcPath,fqFiles[j],qcFolder)
        comm.pools <- paste(comm.pools,"&&",preQC)
        remove.command <- remove.unneeded.files.SE(fqFiles[j])
        comm.pools <- paste(comm.pools,"&&",remove.command)
        move.file.log <- move.paired.files.SE((fqFiles[j]),qcFolder)
        comm.pools <- paste(comm.pools,"&&",move.file.log)

      }
      # Quality control check via Fastqc of result files
      postQC <- postQualityCheck.SE(fastqcPath,fqFiles[j],qcFolder)
      if (is.null(comm.pools)){
        comm.pools <- paste(postQC)
      } else {
        comm.pools <- paste(comm.pools,"&&",postQC)
      }
    }
  }
  if (is.null(comm.pool)){
    comm.pool <- paste(comm.pools,"&")
  } else {
    comm.pool <- paste(comm.pool,comm.pools,"&")
  }
}
store.artificial <- store.artificial.seqs.file(paste(dest.dir,artificial.fq,sep=""),qcFolder)
comm.pool <- paste(comm.pool,"wait","&&", store.artificial)
printOrExecute(comm.pool,Condor)
