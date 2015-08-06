source("GLSeq.Util.R")

# Creates a unique run header
create.text.add <- function(expID,runAttempt) {
  if(is.null(expID) || is.null(runAttempt)) stop("Arguments should not be NULL")
  # the real text.add that will be used as a common ID for all the GLSeq scripts is generated right here:
  expID <- gsub("[[:punct:]]","_",expID)
  expID <- gsub("[[:space:]]","_",expID)
  # The ID + the run attempt is the file and folder naming scheme
  text.add <- paste(expID, runAttempt, sep=".")
  text.add
}

# Destination directory for the processed files (if no connection to DB performed):
create.dest.dir <- function(dest.dir.base,text.add){
  if(is.null(dest.dir.base) || is.null(text.add)) stop("Arguments should not be NULL")
  dest.dir.base <- trailDirCheck(dest.dir.base)
  # Creates the destination directory.  It does so by going to the
  # Destination directory the user input and adding the text.add
  # suffix that is generated for each run
  dest.dir <- paste(dest.dir.base, text.add, "/", sep="")
  dest.dir <- trailDirCheck(dest.dir)
  dest.dir
}

# Destination directory for log /stat files:
create.dest.dir.log <- function(dest.dir,text.add){
  if(is.null(dest.dir) || is.null(text.add)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  # Let's call the status and log files folder .stat
  destDirLog <-  paste(dest.dir, text.add, ".stat/", sep="")
  destDirLog
}

# Instantiates a few of the basic log properties at the start
add.basic.logs <- function(log.file=NULL){
  if(is.null(log.file)) return
  # Just a bunch of logs that are added.
  # This may be expanded/ made customizable in the future to allow for custom loggin
  # These should be sufficient for now.
  add.to.logs(paste("PREPARING DATA:",dataPrepare),log.file)
  add.to.logs(paste("ALIGNING READS:",alignment),log.file)
  add.to.logs(paste("COUNTING ALIGNED READS:",counting),log.file)
  add.to.logs(paste("COLLECTING RESULTS:",resCollect),log.file)
  add.to.logs(paste("PROTOCOL ID:",protID),log.file)
  add.to.logs(paste("ATTRIBUTE FILE LOCATION:",attrPath),log.file)
  add.to.logs("################## END OF RUN OPTIONS ##################",log.file)
  add.to.logs("################## Selected Attribute File Values ##################",log.file)
  add.to.logs(paste("Raw Data Directory:",raw.dir),log.file)
  add.to.logs(paste("Ready Data Directory:",readyData.dir),log.file)
  add.to.logs(paste("Reference FASTA name:",refFASTAname),log.file)
  add.to.logs(paste("Reference GFF name:",refGFFname),log.file)
  add.to.logs(paste("Script Directory:",base.dir),log.file)
  add.to.logs(paste("Destination Directory:",dest.dir),log.file)
  add.to.logs(paste("Alignment Algorithm:",aAlgor),log.file)
  add.to.logs(paste("Counting Algorithm(s):",HTSeq,FeatureCounts,RSEM),log.file)
  add.to.logs("################## END OF ATTRIBUTE OPTIONS ##################",log.file)
}

# Creates a log file
instantiate.logs <- function(dest.dir,text.add,destDirTest=NULL,Condor=FALSE){
  # Creates a log file
  if(is.null(dest.dir) || is.null(text.add)) stop("Arguments should not be NULL")
  log.file <- NULL
  if (!is.null(destDirTest)){
    destDirTest <- trailDirCheck(destDirTest)
    # Creates a run file if there isn't one
    log.file <- paste(destDirTest,text.add,".RunLog.txt",sep="")
    # Overwrites previous file in case run had problems.
    create.log.file <- paste("echo \"RUN LOG FILE\"",">",log.file)
    add.to.logs("############################################################",log.file)
    try(system(create.log.file))
    add.to.logs(paste("Arguments for this run"),log.file)
  } else{
    dest.dir <- trailDirCheck(dest.dir)
    # Creates a log file that has "Run Log File" written
    # on the top via echoing and pushing the stdou to the file
    log.file <- paste(dest.dir,text.add,".RunLog.txt",sep="")
    create.log.file <- paste("echo \"RUN LOG FILE\"",">",log.file)
    printOrExecute(create.log.file,Condor)
  }
  # Handle each case individually in case we end up reworking these more
  try(add.basic.logs(log.file))
  log.file
}

# Creates the run directory (destination directory)
create.run.directory <- function(dest.dir,Condor=FALSE) {
  if(is.null(dest.dir)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  # Creates a new directory with the mkdir command
  destDir.create <- paste("mkdir ", dest.dir, sep="")
  printOrExecute(destDir.create,Condor)
  # Returns command for unit tests
  #destDir.create
}

# Creates a new log directory
create.log.directory <- function(destDirLog,Condor=FALSE){
  if(is.null(destDirLog)) stop("Arguments should not be NULL")
  destDirLog <- trailDirCheck(destDirLog)
  # Creates a new directory with the mkdir command
  destDirLog.create <- paste("mkdir", destDirLog)
  printOrExecute(destDirLog.create,Condor)
  #destDirLog.create
}

# Creates a new HTSeq folder
create.HtSeq.folder <- function (dest.dir, text.add,Condor=FALSE){
  if(is.null(dest.dir) || is.null(text.add)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  # Destination directory for HTSeq
  destDirHTSeqCount <-  paste(dest.dir, text.add, ".HTSeq.Counting/", sep="")
  # Makes a new directory
  destDirHTSeqCount.create <- paste("mkdir ", destDirHTSeqCount, sep="")
  printOrExecute(destDirHTSeqCount.create,Condor)
  destDirHTSeqCount
}

# Creates the folder for FeatureCounts
create.FeatureCounts.folder <- function (dest.dir, text.add,Condor=FALSE){
  if(is.null(dest.dir) || is.null(text.add)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  # Destination directory for FeatureCounts
  destDirFeatureCountsCount <-  paste(dest.dir, text.add, ".FeatureCounts.Counting/", sep="")
  # Just a mkdir command
  destDirFeatureCountsCount.create <- paste("mkdir ", destDirFeatureCountsCount, sep="")
  printOrExecute(destDirFeatureCountsCount.create,Condor)
  destDirFeatureCountsCount
}

# Creates the folder for RSEM
create.RSEM.folder <- function (dest.dir, text.add,Condor=FALSE){
  if(is.null(dest.dir) || is.null(text.add)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  # Creates a new RSEM counting folder
  destDirRSEMCount <-  paste(dest.dir, text.add, ".RSEM.Counting/", sep="")
  # Just a mkdir command
  destDirRSEMCount.create <- paste("mkdir ", destDirRSEMCount, sep="")
  printOrExecute(destDirRSEMCount.create,Condor)
  destDirRSEMCount
}

# Creates the folder for cufflinks results
create.Cufflinks.folder <- function (dest.dir, text.add,Condor=FALSE){
  if(is.null(dest.dir) || is.null(text.add)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  # Creates a new cufflinks results folder
  destDirCufflinksCount <-  paste(dest.dir, text.add, ".Cufflinks.Counting/", sep="")
  # Just making a directory
  destDirCufflinksCount.create <- paste("mkdir ", destDirCufflinksCount, sep="")
  printOrExecute(destDirCufflinksCount.create,Condor)
  destDirCufflinksCount
}

# Create the collection folder normally
create.Collect.folder <- function (dest.dir, text.add,Condor=FALSE){
  if(is.null(dest.dir) || is.null(text.add)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  # Creates a collection folder in the new destination directory
  collectDir <- paste(dest.dir,text.add,".Collect/",sep="")
  # Just a shell mkdir command to make the new directory
  collectDir.create <- paste("mkdir",collectDir)
  printOrExecute(collectDir.create,Condor)
  collectDir
}


# Copy attribute file into the destination directory for record keeping
copy.attribute.file.to.dest <- function(attrPath,dest.dir,Condor=FALSE){
  if(is.null(dest.dir) || is.null(attrPath)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  # Copies the attribute file into the destination directory
  copyAttributeFile <- paste("cp",attrPath,dest.dir)
  printOrExecute(copyAttributeFile,Condor)
  # Returns command for unit tests
  copyAttributeFile
}

# Name some of the log files
create.run.logs <- function(destDirLog,text.add) {
  ###############
  # Log-tail command addition
  ###############
  # Just some naming that we can easily change from here
  runLogFile1 <- paste(destDirLog, text.add, ".runLog1.txt", sep="")
  runLogFile2 <- paste(destDirLog, text.add, ".runLog2.txt", sep="")
  runErrFile1 <- paste(destDirLog, text.add, ".runErr1.txt", sep="")
  runErrFile2 <- paste(destDirLog, text.add, ".runErr2.txt", sep="")
}
############################
# Preparing data file names (ready for expression computation)
############################
# reusable function for assembly split fastq file names (see below)
fqfiles.table.pe.assemble <- function(fqfiles.base) {
  if (is.null(fqfiles.base)) stop("Arguments should not be NULL")
  fqfiles.table <- NULL
  # Add a fq file extension to all of the files
  for (i in 1: length(fqfiles.base)) {
    # This treats all files like .fq file
    left.fq <- paste(fqfiles.base[i],".1.fq", sep="")
    right.fq <- paste(fqfiles.base[i],".2.fq", sep="")
    # Binds the two paired files into a row of a data frame
    fqfiles.table <- rbind(fqfiles.table, c(left.fq, right.fq)) }
  fqfiles.table
}


# Rename processed files
rename.preprocessed.files <- function(data.dir,readyData.dir,Condor=FALSE){
  # Rename .fastq files to .fq files
  if (is.null(dest.dir) || is.null(readyData.dir)) stop("Arguments should not be NULL")
  # Grab the files in the data directory
  # As we won't always be able to look at the files already transfered into the data directory
  setwd(readyData.dir)
  active.files <- dir(readyData.dir)
  setwd(data.dir)
  # Stock starter command
  rename.command <- paste("date")
  for (num in 1:length(active.files)){
    # This reassigns the active file to have the full file path normally
    active.files[num] <- paste(dest.dir,active.files[num],sep="")
    # Only create a rename command for files that end in .fastq
    if (grepl(".fastq$",active.files[num])){
      # Take off the fastq ending
      to.rename.base <- substr(active.files[num],1,nchar(active.files[num])-5)
      # Add the fq ending
      move.file <- paste(to.rename.base,"fq",sep="")
      # Execute that command, simply a move command
      rename.command <- paste(rename.command,";","mv",active.files[num],move.file)
    }
  }
  printOrExecute(rename.command,Condor)
  # Returns command for unit tests
  rename.command
}

# Copying the ready-to-process fastq files to the destination directory:
copy.preprocessed.files <- function(readyData.dir,dest.dir,Condor=FALSE) {
  if (is.null(readyData.dir) || is.null(dest.dir)) stop("Arguments should not be NULL")
  readyData.dir <- trailDirCheck(readyData.dir)
  dest.dir <- trailDirCheck(dest.dir)
  # Creates the copy command that will pull files from the data directory
  # This will specifically take all of the .fq files
  copyFq <- paste(readyData.dir,"*.fq",sep="")
  readyDataCopy <- paste("cp",copyFq,dest.dir)
  # Looks for files in the ready data directory
  active.files <- dir(readyData.dir)
  # If fastq files exist in the dataset
  fastq.files <- FALSE
  # Goes through all the files checking their file extensions
  for (num in 1:length(active.files)){
    if (grepl(".fastq$",active.files[num])){
      # This will cause us to copy all of the .fastq files
      copy.fastq <- paste(readyData.dir,"*.fastq",sep="")
      # Copies the files to the destination directory
      copy.fastq <- paste("cp",copy.fastq,dest.dir)
      readyDataCopy <- paste(readyDataCopy,";",copy.fastq)
      # Notes we have fastq files that need to be renamed.
      fastq.files <- TRUE
    }
  }
  printOrExecute(readyDataCopy,Condor)
  # If there are fastq files we will need to rename them
  if (fastq.files){
    rename.preprocessed.files(dest.dir,Condor)
  }
  # Returns command for unit tests
  readyDataCopy
}

# Converts the input list of files into a table that can then be understood by
# The downstream program to either pair files correctly or to
convert.file.list.to.table <- function(paired.end,readyData.dir) {
  if (is.null(paired.end) || is.null(dest.dir)) stop("Arguments should not be NULL")
  # Paired End
  if (paired.end) {
    fqfiles <- dir(readyData.dir)
    fqfiles <- fqfiles[grep(".fq$|.fastq$", fqfiles)]
    for (num in 1:length(fqfiles)){
      # Remove 8 if a fastq file
      if (grepl(".fastq$",fqfiles[num])){
        fqfiles[num] <- substr(fqfiles.base[num],1,nchar(fqfiles[num]-8))
      # Remove 5 if a fq file
      } else{
        fqfiles[num] <- substr(fqfiles[num], 1,nchar(fqfiles[num]) - 5)
      }
    }
    # Only take the unique ones (This should reduce the number in half)
    fqfiles.base <- unique(fqfiles)
    # Table of pairs of FASTQ file names (1 pair per row):
    # We assemble left and right file names explicitly here because
    # We don't want to rely on file sequence in the directory and get wrong results if there are occasional unpaired files etc.
    # Assemble function for paired-end libraries:
    fqfiles.table <- fqfiles.table.pe.assemble(fqfiles.base)
  }
  # Single End
  if (!(paired.end)) {
    files <- dir(readyData.dir)
    files <- files[grep(".fq$|.fastq$", files)]
    # Fixing the fastq files so they will match when it matters
    for (num in 1:length(files)){
      if (grepl(".fastq$",files[num])){
        # Remove astq
        # Adds a q back to make it a fastq file
        files[num] <- substr(files[num],1,nchar(files[num]-4))
        files[num]<- paste(files[num],"q",sep="")
      }
    }
    fqfiles.table <- cbind(files)
  }
  fqfiles.table
}

# If the data needs to be processed all the way
# from the compressed archives in the raw directory
find.files.for.dataprep <- function(raw.dir,unzipped){
  if (is.null(raw.dir) || is.null(unzipped)) stop("Arguments should not be NULL")
  ##############
  # Paired End
  ##############
  if (paired.end) {
    ########################################
    ############ NO UPDATE #################
    ########################################
    # Raw Directory of files supplied
    gzfiles <- dir(raw.dir)
    fqfiles <- dir(raw.dir)
    if (!unzipped) {
      gzfiles <- gzfiles[grep(".gz", gzfiles)]
      if(!presplit) fqfiles.base <- substr(gzfiles, 1,nchar(gzfiles) - 3)
      if (presplit) fqfiles.base <- substr(gzfiles, 1,nchar(gzfiles) - 3)
      fqfiles.table <- fqfiles.table.pe.assemble(fqfiles.base)
    }
    if (unzipped) {
      fqfiles <- fqfiles[grep(".fq|.fastq",fqfiles)]
      if(!presplit) fqfiles.base <- substr(fqfiles, 1,nchar(fqfiles) - 0)
      if (presplit) fqfiles.base <- substr(fqfiles, 1,nchar(fqfiles) - 0)
      fqfiles.table <- fqfiles.table.pe.assemble(fqfiles.base)
    }
    fqfiles.table
  }
  ##############
  # Single End
  ##############
  if (!(paired.end))  {
    gzfiles <- dir(raw.dir)
    fqfiles <- dir(raw.dir)
    if (!unzipped){
      gzfiles <- dir(raw.dir)
      gzfiles <- gzfiles[grep(".gz", gzfiles)] # just in case raw.dir has somethig else besides .gz files
      fqfiles <- substr(gzfiles, 1, nchar(gzfiles)-3)
      fqfiles.table <- cbind(NULL, fqfiles)
    }
    if (unzipped){
      fqfiles <- dir(raw.dir)
      fqfiles <- fqfiles[grep(".fq|.fastq", fqfiles)] # just in case raw.dir has somethig else besides .fq or .fastq files
      fqfiles <- substr(fqfiles, 1, nchar(fqfiles)-0)
      fqfiles.table <- cbind(NULL, fqfiles)
    }
    fqfiles.table
  }
}

# Chunks the data into streams of roughly equal size
prepare.chunk.function <- function(fqfiles.table,nStreams){
  if (is.null(fqfiles.table) || is.null(nStreams)) stop("Arguments should not be NULL")
  nStreams <- check.nStreams(fqfiles.table,nStreams)
  chunk <- function(x, n) split(x, sort(rank(x) %% n)) # commonly known solution to divide data equally
  rangelist <- chunk(1:nrow(fqfiles.table), nStreams)
  # ranges are in rangelist[1:nSreams]
  for (ii in 1:nStreams) assign(paste("range",ii,sep=""),rangelist[[ii]])
  rangelist
}

# Checks to make sure that the streams are correctly set for the amount of data we have
check.nStreams <- function(fqfiles.table,nStreams){
  # quiet and nice solution but there will de descrepancy between nStreams in the attribute file and the actual nStreams
  if (nStreams > nrow(fqfiles.table)) nStreams <- nrow(fqfiles.table)
  nStreams
}

# Runs the alignment process
start.alignment.process <- function(base.dir,rangelist,nStreams,Condor=FALSE){
  setwd(base.dir)
  this.resName <- NULL
  source("GLSeq.Alignment.R")
  # Returns the big set of commands which is the comm stack
  comm.stack.pool
}

# Runs the counting process if no alignment was run preior
start.counting.process <- function(countable.sams.dir,dest.dir,base.dir,Condor){
  # Heading of the command
  comm.stack.pool <- "date"
  # Takes all the countable.sam files in the folder
  countable.sams <- dir(countable.sams.dir)[grep("countable.sam", dir(countable.sams.dir))]
  for (i in 1:length(countable.sams)){
    countable.sam <- countable.sams[i]
    # Copies
    copy.comm <- paste("cp",paste(countable.sams.dir,countable.sam,sep="/"),dest.dir)
    printOrExecute(copy.comm)
    setwd(base.dir)
    source("GLSeq.Counting.R")
  }
  comm.stack.pool <- paste(comm.stack.pool,"&&",count.comm,"&")

  # RETURNS THE COMM STACK
  comm.stack.pool
}

# When RSEM finishes it needs to do some file movements to make things look nicer
RSEM.finish <- function(comm.stack.pool,destDirRSEMCount,dest.dir) {
  if (is.null(comm.stack.pool) || is.null(destDirRSEMCount) || is.null(dest.dir)) stop("Arguments should not be NULL")
    destDirRSEMCount <- trailDirCheck(destDirRSEMCount)
    dest.dir <- trailDirCheck(dest.dir)
    comm.stack.pool <- paste(comm.stack.pool,"&&","mv",paste(dest.dir,"*.RSEM.counts.*",sep=""),destDirRSEMCount)
    comm.stack.pool <- paste(comm.stack.pool,"&&","mv",paste(dest.dir,"*.index.*",sep=""),destDirRSEMCount)
    # RETURNS THE COMM STACK
    comm.stack.pool
}

# Saves the run data at the end of the run
save.run.data <- function(base.dir,text.add) {
  base.dir <- trailDirCheck(base.dir)
  # Creates an R data file
  currentRun.dataFile <- paste(base.dir,"GLSeq.vars.", text.add, ".rda", sep="")
  if (currentRun.dataFile %in% dir(base.dir)) {
    cat(timestamp(), " Attribute file for this run ID - ", " \n", currentRun.dataFile, " - is already saved; keeping the original file ", "\n")
    currentRun.dataFile.base <- substr(currentRun.dataFile, 1, nchar(currentRun.dataFile) -4)
    add.num <- 1 }
  # saving a version of the attribute file with numeric addon:nStreams
  while(currentRun.dataFile %in% dir(base.dir)) {
    currentRun.dataFile <- paste(currentRun.dataFile.base, add.num, "rda", sep=".")
    add.num <- add.num+1
  }
  if (!(currentRun.dataFile %in% dir(base.dir))) save.image(file=currentRun.dataFile)
  # Returns file name
}

# Runs teh data preparation script
run.data.prep <- function(destDirLog,text.add,attrPath,dest.dir,base.dir,Condor=FALSE){
  if (is.null(destDirLog) || is.null(text.add) || is.null(attrPath) || is.null(dest.dir)) stop("Arguments should not be NULL")
  destDirLog <- trailDirCheck(destDirLog)
  dest.dir <- trailDirCheck(dest.dir)
  Sys.sleep(10)
  dataPrepLog <- paste(destDirLog, text.add, ".DataPrepLog.txt", sep="")
  dataPrepErr <- paste(destDirLog, text.add, ".DataPrepErr.txt", sep="")
  dataPrep <- paste("Rscript ", paste(base.dir,GLSeq.dataprep.R,sep=""), text.add," ",dest.dir," ",attrPath," 1>> ", dataPrepLog, " 2>> ", dataPrepErr, sep="")
  printOrExecute(dataPrep,Condor)
}

# Executes the communication stack
execute.comm.stack <- function (comm.stack.pool, Condor=FALSE) {
  if (is.null(comm.stack.pool)) stop("Arguments should not be NULL")
  stack.start.time <- proc.time()
  printOrExecute(comm.stack.pool,Condor)
}

# Creates the collection directories
previous.run.directories <- function(previous.dir,previous.run.name){
  if (is.null(previous.dir) || is.null(previous.run.name)) stop("Arguments should not be NULL")
  dest.dir <- previous.dir
  dest.dir <- trailDirCheck(dest.dir)
  dest.dir
}
previous.run.FeatureCounts <- function(previous.dir,previous.run.name){
  if (is.null(previous.dir) || is.null(previous.run.name)) stop("Arguments should not be NULL")
  dest.dir <- previous.run.directories(previous.dir,previous.run.name)
  dest.dir <- trailDirCheck(dest.dir)
  destDirFeatureCountsCount <-  paste(dest.dir,previous.run.name,".FeatureCounts.counting/", sep="")
  destDirFeatureCountsCount
}
previous.run.HTSeq <- function(previous.dir,previous.run.name){
  if (is.null(previous.dir) || is.null(previous.run.name)) stop("Arguments should not be NULL")
  dest.dir <- previous.run.directories(previous.dir,previous.run.name)
  dest.dir <- trailDirCheck(dest.dir)
  destDirHTSeqCount <-  paste(dest.dir,previous.run.name,".HTSeq.counting/", sep="")
  destDirHTSeqCount
}
previous.run.RSEM <- function(previous.dir,previous.run.name){
  if (is.null(previous.dir) || is.null(previous.run.name)) stop("Arguments should not be NULL")
  dest.dir <- previous.run.directories(previous.dir,previous.run.name)
  dest.dir <- trailDirCheck(dest.dir)
  destDirRSEMCount <-  paste(dest.dir,previous.run.name,".RSEM.counting/", sep="")
  destDirRSEMCount
}
collect.counted.results <- function() {
  setwd(base.dir)
  source("GLSeq.ResultsCollect.R")
}