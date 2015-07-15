source("../Main/GLSeq.Util.R")

###########################
# loading variables for the current run (from the current directory) 
###########################
#
load.dataFile <- function(text.add){
  if (is.null(text.add)) stop("Arguments should not be NULL")
  currentRun.dataFile <- paste("GLSeq.vars.", text.add, ".rda", sep="")
  try(load(currentRun.dataFile)) # loads all parameters for the run with the runID indicated in the text.add object
  currentRun.dataFile
}
#
###########################
# assembling Trimmomatic system commands: 
###########################
trimAssemble.PE <- function(fqFile, trimPath, qScore, headcrop=12, artificial.fq,trimMin) {
  if (is.null(fqFile) || is.null(trimPath) || is.null(qScore) || is.null(artificial.fq) || is.null(trimMin)){
    stop("Arguments should not be NULL")
  }
  if (substr(trimPath,nchar(trimPath),nchar(trimPath))=="/") stop("Invalid path to Trimmomatic JAR file")
  fqFile <- get.file.base(fqFile)
  first.read.file <- first.read.name(fqFile)
  second.read.file <- second.read.name(fqFile)
  
  dashqScore <- paste("-", qScore, sep="")
  #
  logName <- paste(first.read.file,"pairedtrim.log", sep=".")
  #
  pairedTrimmed.1 <- paste("p", first.read.file, sep=".")
  unpairedTrimmed.1 <- paste("u",  first.read.file, sep=".")
  #
  pairedTrimmed.2 <- paste("p",  second.read.file, sep=".")
  unpairedTrimmed.2 <- paste("u", second.read.file, sep=".") 
  #
  trimParam <- paste("ILLUMINACLIP:", artificial.fq, ":2:30:10", " HEADCROP:", headcrop, " SLIDINGWINDOW:3:30 MINLEN:", trimMin, sep="")
  tcomm <- paste("java -jar", trimPath, "PE -threads 20", dashqScore, "-trimlog", logName, first.read.file, second.read.file, pairedTrimmed.1, unpairedTrimmed.1, pairedTrimmed.2, unpairedTrimmed.2, trimParam)
  # Return
  tcomm
}

trimAssemble.SE <- function(fqFile, trimPath, qScore, headcrop=12, artificial.fq, trimMin) {
  if (is.null(fqFile) || is.null(trimPath) || is.null(qScore) || is.null(artificial.fq) || is.null(trimMin)){
    stop("Arguments should not be NULL")
  }
  if (substr(trimPath,nchar(trimPath),nchar(trimPath))=="/") stop("Invalid path to Trimmomatic JAR file")
  dashqScore <- paste("-", qScore, sep="")
  logName <- paste(fqFile, "pairedtrim.log", sep=".")
  pairedTrimmed.1 <- paste("p", fqFile, "fq", sep=".") # '.fq' is added here for single end libraries (it is added at the splitting stage for the paired-end)
  trimParam <- paste("ILLUMINACLIP:", artificial.fq, ":2:30:10", " HEADCROP:", headcrop, " SLIDINGWINDOW:3:30 MINLEN:", trimMin, sep="")
  tcomm <- paste("java -jar", trimPath, "SE -threads 20", dashqScore, "-trimlog", logName, fqFile,  pairedTrimmed.1, trimParam)
  tcomm
}

get.files.unzipped <- function(raw.dir){
  if (is.null(raw.dir)) stop("Arguments should not be NULL")
  raw.dir <- trailDirCheck(raw.dir)
  if (length(unique(raw.dir)) == 1) {
    # Now takes either .fq or .fastq
    fqFiles.unzip <- dir(raw.dir[1])[grep(".fq$|.fastq$", dir(raw.dir[1]))]
  }
  fqFiles.unzip
}

get.files.zipped <- function(raw.dir){
  if (is.null(raw.dir)) stop("Arguments should not be NULL")
  raw.dir <- trailDirCheck(raw.dir)
  # Gets the zipped files with the compressed extension .gz
  if (length(unique(raw.dir)) == 1) {
    fqFiles.zip <- dir(raw.dir[1])[grep(".fq.gz$|.fastq.gz$", dir(raw.dir[1]))]
  }
  fqFiles.zip
}


prepare.unzipped.file.names <- function(fqFiles.unzip){
  if (is.null(fqFiles.unzip))  stop("Arguments should not be NULL")
  # Assumes the file is in the format of {name}_#.fq (Presplit) OR {name}_fq
  # Thus, the file example.1.fq along with its pair of example.2.fq
  # would become example.1.fq and example.2.fq
  # and an unsplit file example.fq would become example.fq
  fqFiles <- substr(fqFiles.unzip, 1, nchar(fqFiles.unzip)-0)
  fqFiles
}

prepare.zipped.file.names <- function(fqFiles.zip){
  if (is.null(fqFiles.zip))  stop("Arguments should not be NULL")
  # Assumes file is in the format of {name}_#.fq.gz (Presplit) OR {name}.fq.gz (Unsplit)
  # Thus, the file example.1.fq.gz along with its pair of example.2.fq.gz
  # would become example.1.fq and example.2.fq
  # and an unsplit file example.fq.gz would become example.fq
  fqFiles <- substr(fqFiles.zip, 1, nchar(fqFiles.zip)-3)
  fqFiles
}

check.presplit <- function(paired.end,presplit){
  if (is.null(paired.end) || is.null(presplit))  stop("Arguments should not be NULL")
  if (!is.logical(paired.end) || !is.logical(presplit)) stop("Arguments must be LOGICALS")
  if (!paired.end) presplit <- FALSE
  presplit
}

check.nStreamsDataPrep <- function(fqFiles,nStreamsDataPrep){
  if (is.null(fqFiles) || is.null(nStreamsDataPrep))  stop("Arguments should not be NULL")
  if (nStreamsDataPrep > nrow(fqFiles)) nStreamsDataPrep <- nrow(fqFiles)
  nStreamsDataPrep
}

chunk.data.files.presplit <- function(fqFiles,nStreamsDataPrep){
  if (is.null(fqFiles) || is.null(nStreamsDataPrep))  stop("Arguments should not be NULL")
  nStreamsDataPrep <- check.nStreamsDataPrep(fqFiles,nStreamsDataPrep)
  chunk <- function(x, n) split(x, sort(rank(x) %% n)) # commonly known solution to divide data equally
  rangelist.Dataprep <- chunk(2*(1:(nrow(fqFiles)/2)), nStreamsDataPrep)
  rangelist.Dataprep
}

chunk.data.files.unsplit <- function(fqFiles,nStreamsDataPrep){
  if (is.null(fqFiles) || is.null(nStreamsDataPrep))  stop("Arguments should not be NULL")
  nStreamsDataPrep <- check.nStreamsDataPrep(fqFiles,nStreamsDataPrep)
  chunk <- function(x, n) split(x, sort(rank(x) %% n)) # commonly known solution to divide data equally
  rangelist.Dataprep <- chunk(1:nrow(fqFiles), nStreamsDataPrep)
  rangelist.Dataprep
}

check.ifFiles <- function(fqFiles.zip,fqFiles.unzip){
  if (is.null(fqFiles.zip) && is.null(fqFiles.unzip)) stop("Please check the list of raw FASTQ files and the contents of the raw directory \n Please remember that files must be in .fq or .fq.gz format \n")
  NULL
}

copy.artificial.fq <- function(base.dir,artificial.fq,dest.dir){
  if (is.null(base.dir) || is.null(artificial.fq) || is.null(dest.dir)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  base.dir <- trailDirCheck(base.dir)
  copy <- paste("cp",paste(base.dir,artificial.fq,sep=""),dest.dir)
  try(system(copy))
  copy
}

get.file.base <- function(fqFile){
  if (grepl(".gz$",fqFile)){
    fqFile <- substr(fqFile, 1, nchar(fqFile)-3)
  }
  if (grepl(".fastq$",fqFile)){
    fqFile <- substr(fqFile, 1, nchar(fqFile)-6)
  }
  if (grepl(".fq$",fqFile)){
    fqFile <- substr(fqFile, 1, nchar(fqFile)-3)
  }
  if (grepl(".1$",fqFile)){
    fqFile <- substr(fqFile, 1, nchar(fqFile)-2)
  }
  if (grepl(".2$",fqFile)){
    fqFile <- substr(fqFile, 1, nchar(fqFile)-2)
  }
  fqFile
}

copy.file<- function(raw.dir,fqFile,dest.dir){
  if (is.null(raw.dir) || is.null(fqFile) || is.null(dest.dir)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  raw.dir <- trailDirCheck(raw.dir)
  #
  copy.comm <- paste("cp",paste(raw.dir,fqFile,sep=""),dest.dir)
  copy.comm
}

unzip.gz.files <- function(fqFile.zip){
  if (is.null(fqFile.zip)) stop("Arguments should not be NULL")
  #
  gunzip.comm <- paste("gunzip",fqFile.zip)
  gunzip.comm
}

first.read.name <- function(fqFile){
  if (is.null(fqFile)) stop("Arguments should not be NULL")
  fqFile <- get.file.base(fqFile)
  first.read.filename <- paste(fqFile, ".1.fq", sep="")
  first.read.filename
}

second.read.name <- function(fqFile){
  if (is.null(fqFile)) stop("Arguments should not be NULL")
  fqFile <- get.file.base(fqFile)
  second.read.filename <- paste(fqFile, ".2.fq", sep="")
  second.read.filename
}

split.unsplit.files.PE <- function(dest.dir,fqFile){
  if (is.null(dest.dir) || is.null(fqFile)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  fqFile.base <- get.file.base(fqFile)
  first.read.filename <- first.read.name(fqFile.base)
  second.read.filename <- second.read.name(fqFile.base)
  # (The ruby code incorporated into the line above was adapted from SeqAnsweres forum (seqanswers.com))
  split.comm <- paste("cat ", dest.dir, fqFile, " | ruby -ne 'BEGIN{@i=0} ; @i+=1; puts $_  if @i.to_s =~ /[1234]/; @i = 0 if @i == 8' > ", first.read.filename, " && cat ",   dest.dir, fqFile, " | ruby -ne 'BEGIN{@i=0} ; @i+=1; puts $_  if @i.to_s =~ /[5678]/; @i = 0 if @i == 8' > ", second.read.filename, sep="")
  split.comm
}

left.dirty.name <- function(fqFile){
  if (is.null(fqFile)) stop("Arguments should not be NULL")
  fqFile <- get.file.base(fqFile)
  leftDirtyFname <-  paste("dirty.", fqFile, ".1.fq", sep="")
  leftDirtyFname
}

right.dirty.name <- function(fqFile){
  if (is.null(fqFile)) stop("Arguments should not be NULL")
  fqFile <- get.file.base(fqFile)
  rightDirtyFname <-  paste("dirty.", fqFile, ".2.fq", sep="")
  rightDirtyFname
}

#
#
# Paired Ended Functions
#
#
#

file.shuffle.PE <- function(fqFile){
  if (is.null(fqFile)) stop("Arguments should not be NULL")
  fqFile <- get.file.base(fqFile)
  first.read.filename <- first.read.name(fqFile)
  second.read.filename <- second.read.name(fqFile)
  leftDirtyFname <- left.dirty.name(fqFile)
  rightDirtyFname <- right.dirty.name(fqFile)
  pairedTrimmed.1 <-  paste("p.", fqFile, ".1.fq", sep="")
  pairedTrimmed.2 <- paste("p.", fqFile, ".2.fq", sep="")
  #
  fileShuffle <- paste("mv", first.read.filename, leftDirtyFname, "&&", "mv", second.read.filename, rightDirtyFname, "&&", "mv", pairedTrimmed.1, first.read.filename, "&&", "mv", pairedTrimmed.2, second.read.filename)
  #
  #
  fileShuffle
}

preQualityCheck.PE <- function(fastqcPath, fqFile){
  if (is.null(fastqcPath) || is.null(fqFile)) stop("Arguments should not be NULL")
  fqFile <- get.file.base(fqFile)
  leftDirtyFname <- left.dirty.name(fqFile)
  rightDirtyFname <- right.dirty.name(fqFile)
  preQC <- paste(fastqcPath, leftDirtyFname, rightDirtyFname)
  preQC
}

postQualityCheck.PE <- function(fastqcPath,fqFile){
  fqFile <- get.file.base(fqFile)
  if (is.null(fastqcPath) || is.null(fqFile)) stop("Arguments should not be NULL")
  first.read.filename <- first.read.name(fqFile)
  second.read.filename <- second.read.name(fqFile)
  postQC <- paste(fastqcPath, first.read.filename, second.read.filename)
  postQC
}