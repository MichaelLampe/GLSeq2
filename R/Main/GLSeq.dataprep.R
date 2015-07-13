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
attrPath <- as.character(args[2])
source("GLSeq.Util.R")
source(attrPath)
#
###########################
# loading variables for the current run (from the current directory) 
###########################
#
currentRun.dataFile <- paste("GLSeq.vars.", text.add, ".rda", sep="")
load(currentRun.dataFile) # loads all parameters for the run with the runID indicated in the text.add object
#
################################################################################
##################### PROCESS COMPRESSED FASTQ FILES ###########################
################################################################################
# Routine for compressed (".gz") fastq files, paired-end sequencing and 
# concatenated (1-st and 2-nd reads pooled) FASTQ files,  
# i.e. current situation in GLBRC/JGI for both E.coli and yeast data (April 2013);
# with appearence of new data types, respective data preparation blocks will 
# be added as needed
################################################################################
# 
###########################
# assembling Trimmomatic system commands: 
###########################
#
if (paired.end) {
  trimAssemble <- function(leftDirtyName, rightDirtyName, trimPath, qScore, headcrop=12, artifactsFile) {
    #
    dashqScore <- paste("-", qScore, sep="")
    #
    logName <- paste(leftDirtyName, "pairedtrim.log", sep=".")
    #
    pairedTrimmed.1 <- paste("p", leftDirtyName, sep=".")
    unpairedTrimmed.1 <- paste("u",  leftDirtyName, sep=".")
    #
    pairedTrimmed.2 <- paste("p",  rightDirtyName, sep=".")
    unpairedTrimmed.2 <- paste("u", rightDirtyName, sep=".") 
    #
    trimParam <- paste("ILLUMINACLIP:", artifactsFile, ":2:30:10", " HEADCROP:", headcrop, " SLIDINGWINDOW:3:30 MINLEN:", trimMin, sep="")
    tcomm <- paste("java -jar", trimPath, "PE -threads 20", dashqScore, "-trimlog", logName, leftDirtyName, rightDirtyName, pairedTrimmed.1, unpairedTrimmed.1, pairedTrimmed.2, unpairedTrimmed.2, trimParam)
    # Return
    tcomm
  }
}
#
###########################
# With single-end libraries, Trimmomatic commands are shorter.
# Let's keep 'p' prefix for the trimmed .fq files for consistency
###########################
#
if (!(paired.end)) {
  trimAssemble <- function(leftDirtyName, trimPath, qScore, headcrop=12, artifactsFile) {
    dashqScore <- paste("-", qScore, sep="")
    logName <- paste(leftDirtyName, "pairedtrim.log", sep=".")
    pairedTrimmed.1 <- paste("p", leftDirtyName, "fq", sep=".") # '.fq' is added here for single end libraries (it is added at the splitting stage for the paired-end)
    trimParam <- paste("ILLUMINACLIP:", artifactsFile, ":2:30:10", " HEADCROP:", headcrop, " SLIDINGWINDOW:3:30 MINLEN:", trimMin, sep="")
    tcomm <- paste("java -jar", trimPath, "SE -threads 20", dashqScore, "-trimlog", logName, leftDirtyName,  pairedTrimmed.1, trimParam)
    tcomm
  }
}
###########################
setwd(dest.dir)
###########################
# Files indicating readiness of the libraries (i.e. split FQ files):
###########################
files2watch.dataprep <- NULL
fqFiles.zip <- NULL
fqFiles.unzip <- NULL
###########################
#
#
#
################################################################################
##################### GETTING COMPRESSED FQ FILE LIST ##########################
################################################################################
#
###########################
# 1) Using supplied list of raw files (planned to be the standard way): 
###########################
#
if (unzipped){
  if (!(is.null(libList))) {
    fqFiles.unzip <- libList
  }
}
#
if (!unzipped){
  if (!(is.null(libList))) fqFiles.zip <- libList
}
#
###########################
# 2) All the raw files are in the same directory 
# and the list of files is not explicitly supplied (the situation we had before introducing a route of communication with CBDB): 
###########################
#
if (unzipped){  
  if (length(unique(raw.dir)) == 1){
    # Now takes either .fq or .fastq
    fqFiles.unzip <- dir(raw.dir[1])[grep(".fq|fastq", dir(raw.dir[1]))]
  }
}
#
#
if (!unzipped){
  if (length(unique(raw.dir)) == 1 & is.null(libList))  fqFiles.zip <- dir(raw.dir[1])[grep("fq.gz|fastq.gz", dir(raw.dir[1]))]
}
#
############################
# If fqFiles.zip and fqFiles.unzip vectors are still empty for some reason:
# Kill program with error
###########################
#
if (is.null(fqFiles.zip) && is.null(fqFiles.unzip)) stop("Please check the list of raw FASTQ files and the contents of the raw directory \n Please remember that files must be in .fq or .fq.gz format \n")
#
###########################
# Names of the respective uncompressed FASTQ files (to be generated): 
# Takes whole file name minus the ''.fq.gz' or '.fq' endings
###########################
#
if (unzipped){
  # Assumes the file is in the format of {name}_#.fq (Presplit) OR {name}_fq
  # Thus, the file example.1.fq along with its pair of example.2.fq
  # would become example.1.fq and example.2.fq
  # and an unsplit file example.fq would become example.fq
  fqFiles <- substr(fqFiles.unzip, 1, nchar(fqFiles.unzip)-0)
}
#
if (!unzipped){
  # Assumes file is in the format of {name}_#.fq.gz (Presplit) OR {name}.fq.gz (Unsplit)
  # Thus, the file example.1.fq.gz along with its pair of example.2.fq.gz
  # would become example.1.fq and example.2.fq
  # and an unsplit file example.fq.gz would become example.fq
  fqFiles <- substr(fqFiles.zip, 1, nchar(fqFiles.zip)-3)
}

#
###########################
# Ranges for data preparation: 
###########################
# Special case to push paired end into !presplit group
# This is to make it so the user doesn't have to care about the presplit command
# When working with SE data
#
if (nStreamsDataPrep > length(fqFiles)) nStreamsDatapPep <- length(fqFiles)
if (!paired.end) presplit <- FALSE
chunk <- function(x, n) split(x, sort(rank(x) %% n)) # commonly known solution to divide data equally
if (presplit) rangelist.Dataprep <- chunk(2*(1:(length(fqFiles)/2)), nStreamsDataPrep)
if (!presplit) rangelist.Dataprep <- chunk(1:length(fqFiles), nStreamsDataPrep)
#
for (zz in 1:nStreamsDataPrep) {
  if (zz==1) comm.pool <- "date"
  if (zz!=1) comm.pool <- paste(comm.pool,"date")
  #
  ############################################################################################################
  ######################################## UNSPLIT PRESPLIT FILES ############################################
  ############################################################################################################
  #
  # These are zipped or unzipped files created
  # From paired-end jobs that have already been
  # Divided into two files
  #
  if (presplit){
    for (j in rangelist.Dataprep[[zz]]) {
      copy.comm <- "date"
      #
      if (j==2) copy.comm <- paste(copy.comm,"&&","cp",paste(base.dir,artificial.fq,sep=""),dest.dir)
      ###########################
      # Copy and Unzip files (If necessary)
      ###########################
      #
      # Gunzip is what does the unzipping, just use cp to copy.
      #
      if (unzipped){
        copy.comm <- paste(copy.comm,"&&","cp",paste(raw.dir,fqFiles.unzip[j-1],sep=""),dest.dir,";","cp",paste(raw.dir,fqFiles.unzip[j],sep=""),dest.dir)
      }
      if(!unzipped){
        copy.comm <- paste(copy.comm,"&&","cp",paste(raw.dir,fqFiles.zip[j-1],sep=""),dest.dir,";","cp",paste(raw.dir,fqFiles.zip[j],sep=""),dest.dir)
        gunzip.comm <- paste("gunzip", fqFiles.zip[j-1],";","gunzip",fqFiles.zip[j])
      }
      fqFile.base <- fqFiles
      first.read.filename <- paste(fqFile.base[j-1], ".fq", sep="")
      second.read.filename <- paste(fqFile.base[j], ".fq", sep="")
      #
      ###########################
      # Trimming Reads
      ###########################
      #
      if (readTrim) {
        leftDirtyFname <-  paste("dirty.", fqFile.base[j-1], ".fq", sep="")
        rightDirtyFname <-  paste("dirty.", fqFile.base[j], ".fq", sep="")
        pairedTrimmed.1 <-  paste("p.", fqFile.base[j-1], ".fq", sep="")
        pairedTrimmed.2 <- paste("p.", fqFile.base[j], ".fq", sep="")
        trimCommand <- trimAssemble(first.read.filename, second.read.filename, trimPath, qScores, trimhead, artificial.fq)
        fileShuffle <- paste("mv", first.read.filename, leftDirtyFname, "&&", "mv", second.read.filename, rightDirtyFname, "&&", "mv", pairedTrimmed.1, first.read.filename, "&&", "mv", pairedTrimmed.2, second.read.filename)
        preQC <- paste(fastqcPath, leftDirtyFname, rightDirtyFname)
        postQC <- paste(fastqcPath, first.read.filename, second.read.filename)
      }
      #
      if (!readTrim){
        preQC <- paste(fastqcPath, first.read.filename, second.read.filename)
      }
      #
      ###########################
      # Construct Command Stack
      ###########################
      # Put these in order
      #
      comm.pool <- paste(comm.pool,"&&",copy.comm)
      if (!unzipped) comm.pool <- paste(comm.pool,"&&",gunzip.comm)
      if (readTrim) comm.pool <- paste(comm.pool,"&&",trimCommand,"&&",fileShuffle)
      comm.pool <- paste(comm.pool,"&&",preQC)
      if(readTrim) comm.pool <- paste(comm.pool,"&&",postQC)
      comm.pool <- paste(comm.pool)
    }
  }
  #
  ############################################################################################################
  #################################### UNSPLIT PAIRED-END FILES ##############################################
  ############################################################################################################
  #
  # These are zipped or unzipped files created
  # From paired-end jobs that will be split into two
  # And then processed accordingly
  #
  if (!presplit && paired.end){
    for (j in rangelist.Dataprep[[zz]]) {
      copy.comm <- "date"
      #
      ####################
      # UNZIPPED
      ####################
      if (unzipped){
        if (j==1){
          copy.comm <- paste(copy.comm,"&&","cp",paste(base.dir,artificial.fq,sep=""),dest.dir)
          copy.comm <- paste(copy.comm,"&&","cp",paste(raw.dir,fqFiles.unzip[j],sep=""),dest.dir)
        }
        if (j!=1) copy.comm <- paste(copy.comm,"&&","cp",paste(raw.dir,fqFiles.unzip[j],sep=""),dest.dir)
      }
      #
      ####################
      # ZIPPED
      ####################
      if (!unzipped){
        if (j==1){
          copy.comm <- paste(copy.comm,"&&","cp",paste(base.dir,artificial.fq,sep=""),dest.dir)
          copy.comm <- paste(copy.comm,"&&","cp",paste(raw.dir,fqFiles.zip[j],sep=""),dest.dir)
        }
        if (j!=1) copy.comm <- paste(copy.comm,"&&","cp",paste(raw.dir,fqFiles.zip[j],sep=""),dest.dir)
        gunzip.comm <- paste("gunzip",fqFiles.zip[j])
        copy.comm <- paste(copy.comm,"&&",gunzip.comm)
      }
      #
      fqFile.base <- fqFiles
      first.read.filename <- paste(fqFile.base[j], ".1.fq", sep="")
      second.read.filename <- paste(fqFile.base[j], ".2.fq", sep="")
      #
      ##################
      # SPLITTING
      ##################
      # (The ruby code incorporated into the line above was adapted from SeqAnsweres forum (seqanswers.com))
      split.comm <- paste("cat ", dest.dir, fqFiles[j], " | ruby -ne 'BEGIN{@i=0} ; @i+=1; puts $_  if @i.to_s =~ /[1234]/; @i = 0 if @i == 8' > ", first.read.filename, " && cat ",   dest.dir, fqFiles[j], " | ruby -ne 'BEGIN{@i=0} ; @i+=1; puts $_  if @i.to_s =~ /[5678]/; @i = 0 if @i == 8' > ", second.read.filename, sep="")
      #
      ##################
      # READ TRIM
      ##################
      if (readTrim) {
        leftDirtyFname <-  paste("dirty.", fqFile.base[j], ".1.fq", sep="")
        rightDirtyFname <-  paste("dirty.", fqFile.base[j], ".2.fq", sep="")
        pairedTrimmed.1 <-  paste("p.", fqFile.base[j], ".1.fq", sep="")
        pairedTrimmed.2 <- paste("p.", fqFile.base[j], ".2.fq", sep="")
        #
        trimCommand <- trimAssemble(first.read.filename, second.read.filename, trimPath, qScores, trimhead, artificial.fq)
        #
        fileShuffle <- paste("mv", first.read.filename, leftDirtyFname, " && ", "mv", second.read.filename, rightDirtyFname, " && ", "mv", pairedTrimmed.1, first.read.filename, " && ", "mv", pairedTrimmed.2, second.read.filename)
        #
        preQC <- paste(fastqcPath, leftDirtyFname, rightDirtyFname)
        #
        postQC <- paste(fastqcPath, first.read.filename, second.read.filename)
        #
        comm.pool <- paste(comm.pool,"&&",copy.comm,"&&",split.comm,"&&",trimCommand,"&&",fileShuffle,"&&",preQC,"&&",postQC)
      }
      #
      ##################
      # NO READ TRIM
      ##################
      if (!readTrim){
        preQC <- paste(fastqcPath, first.read.filename, second.read.filename)
        comm.pool <- paste(comm.pool,"&&",copy.comm,"&&",split.comm,"&&",preQC)
      }  
    }
  }
  ############################################################################################################
  ######################################## SINGLE ENDED FILES ################################################
  ############################################################################################################
  #
  if(!paired.end){
    for (j in rangelist.Dataprep[[zz]]) {
      copy.comm <- "date"
      #
      #
      ####################
      # UNZIPPED
      ####################
      if (unzipped){
        if (j==1) {
          copy.comm <- paste(copy.comm,"&&","cp",paste(base.dir,artificial.fq,sep=""),dest.dir)
          copy.comm <- paste(copy.comm,"&&","cp",paste(raw.dir,fqFiles.unzip[j],sep=""),dest.dir)
        }
        if (j!=1) copy.comm <- paste(copy.comm,"&&","cp",paste(raw.dir,fqFiles.unzip[j],sep=""),dest.dir)
        comm.pool <- paste(comm.pool,"&&",copy.comm)
        fqFile.base <- fqFiles
      }
      #
      ####################
      # ZIPPED
      ####################
      if (!unzipped){
        if (j==1) {
          copy.comm <- paste(copy.comm,"&&","cp",paste(base.dir,artificial.fq,sep=""),dest.dir)
          copy.comm <- paste(copy.comm,"&&","cp",paste(raw.dir,fqFiles.zip[j],sep=""),dest.dir)
        }
        if (j!=1) copy.comm <- paste(copy.comm,"&&","cp",paste(raw.dir,fqFiles.zip[j],sep=""),dest.dir)
        gunzip.comm <- paste("gunzip",fqFiles.zip[j])
        comm.pool <- paste(comm.pool,"&&",copy.comm,"&&",gunzip.comm)
        fqFile.base <- fqFiles
      }
      #
      ####################
      # Trimming
      ####################
      if(readTrim){
        unpaired.fq <- paste(fqFile.base[j], ".fq", sep="")
        SE.dirtyFname <-  paste("dirty.", fqFile.base[j],".fq", sep="")
        SE.trimmedFname <-  paste("p.", unpaired.fq,".fq",sep="")
        #
        trimCommand <- trimAssemble(unpaired.fq, trimPath, qScores, trimhead, artificial.fq)
        #
        fileShuffle <- paste("mv", unpaired.fq, SE.dirtyFname, " && ",  "mv", SE.trimmedFname, unpaired.fq)
        #
        preQC <- paste(fastqcPath, SE.dirtyFname)
        #
        postQC <- paste(fastqcPath, unpaired.fq)
        # moba
        comm.pool <- paste(comm.pool,"&&",trimCommand,"&&",fileShuffle,"&&",preQC,"&&",postQC)
      }
      #
      ####################
      # Not Trimming
      ####################
      if (!readTrim){
        unpaired.fq <- paste(fqFile.base[j], ".fq", sep="") 
        preQC <- paste(fastqcPath, unpaired.fq)
        comm.pool <- paste(comm.pool,"&&",preQC)
      }
    }
    comm.pool <- paste(comm.pool,"&")
  }
}
dataReady.signal <- paste("echo > ", text.add, ".DataReady", sep="") 
comm.pool <- paste(comm.pool,"wait","&&",dataReady.signal)
system(comm.pool)