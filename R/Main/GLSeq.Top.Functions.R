source("GLSeq.Util.R")

create.text.add <- function(expID,runAttempt) {
  if(is.null(expID) || is.null(runAttempt)) stop("Arguments should not be NULL")
  # the real text.add that will be used as a common ID for all the GLSeq scripts is generated right here:
  expID <- gsub("[[:punct:]]","_",expID)
  expID <- gsub("[[:space:]]","_",expID)
  #
  text.add <- paste(expID, runAttempt, sep=".")
  text.add
}
# Destination directory for the processed files (if no connection to DB performed):
create.dest.dir <- function(dest.dir.base,text.add){
  if(is.null(dest.dir.base) || is.null(text.add)) stop("Arguments should not be NULL")
  dest.dir.base <- trailDirCheck(dest.dir.base)
  dest.dir <- paste(dest.dir.base, text.add, "/", sep="")
  dest.dir <- trailDirCheck(dest.dir)
  dest.dir
}
# Destination directory for log /stat files:
create.dest.dir.log <- function(dest.dir,text.add){
  if(is.null(dest.dir) || is.null(text.add)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  destDirLog <-  paste(dest.dir, text.add, ".stat/", sep="")
  destDirLog
}
add.basic.logs <- function(log.file=NULL){
  if(is.null(log.file)) stop("Log file is null")
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
  if (aAlgor == "Cushaw") {
    if (GPU.accel){
      add.to.logs("Alignment Algorithm: Cushaw-GPU",log.file)
    }else{
      add.to.logs(paste("Alignment Algorithm:",aAlgor),log.file)
    }
  } else{
    add.to.logs(paste("Alignment Algorithm:",aAlgor),log.file)
  }
  add.to.logs(paste("Counting Algorithm(s):",HTSeq,FeatureCounts,RSEM),log.file)
  add.to.logs("################## END OF ATTRIBUTE OPTIONS ##################",log.file)
}
instantiate.logs <- function(dest.dir,text.add,destDirTest=NULL){
  # Creates a log file
  if(is.null(dest.dir) || is.null(text.add)) stop("Arguments should not be NULL")
  log.file <- NULL
  if (!is.null(destDirTest)){
    destDirTest <- trailDirCheck(destDirTest)
    log.file <- paste(destDirTest,text.add,".RunLog.txt",sep="")
    # Overwrites previous file in case run had problems.
    create.log.file <- paste("echo \"RUN LOG FILE\"",">",log.file)
    add.to.logs("############################################################",log.file)
    try(system(create.log.file))
    add.to.logs(paste("Arguments for this run"),log.file)
  } else{
    dest.dir <- trailDirCheck(dest.dir)
    log.file <- paste(dest.dir,text.add,".RunLog.txt",sep="")
    create.log.file <- paste("echo \"RUN LOG FILE\"",">",log.file)
    try(system(create.log.file))
  }
  # Handle each case individually in case we end up reworking these more
  try(add.basic.logs(log.file))
  log.file
}
create.run.directory <- function(dest.dir) {
  if(is.null(dest.dir)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  destDir.create <- paste("mkdir ", dest.dir, sep="")
  try(system(destDir.create))
  #
  # Returns command for unit tests
  destDir.create
}
create.log.directory <- function(destDirLog){
  if(is.null(destDirLog)) stop("Arguments should not be NULL")
  destDirLog <- trailDirCheck(destDirLog)
  destDirLog.create <- paste("mkdir", destDirLog)
  try(system(destDirLog.create))
  destDirLog.create
}
create.HtSeq.folder <- function (dest.dir, text.add,log.file=NULL){
  if(is.null(dest.dir) || is.null(text.add)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  destDirHTSeqCount <-  paste(dest.dir, text.add, ".HTSeq.Counting/", sep="")
  destDirHTSeqCount.create <- paste("mkdir ", destDirHTSeqCount, sep="")
  #
  add.to.logs(destDirHTSeqCount.create,log.file)
  try(system(destDirHTSeqCount.create))
  destDirHTSeqCount
}
create.FeatureCounts.folder <- function (dest.dir, text.add,log.file=NULL){
  if(is.null(dest.dir) || is.null(text.add)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  destDirFeatureCountsCount <-  paste(dest.dir, text.add, ".FeatureCounts.Counting/", sep="")
  destDirFeatureCountsCount.create <- paste("mkdir ", destDirFeatureCountsCount, sep="")
  #
  add.to.logs(destDirFeatureCountsCount.create,log.file)
  try(system(destDirFeatureCountsCount.create))
  destDirFeatureCountsCount
}
create.RSEM.folder <- function (dest.dir, text.add,log.file=NULL){
  if(is.null(dest.dir) || is.null(text.add)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  destDirRSEMCount <-  paste(dest.dir, text.add, ".RSEM.Counting/", sep="")
  destDirRSEMCount.create <- paste("mkdir ", destDirRSEMCount, sep="")
  #
  add.to.logs(destDirRSEMCount.create,log.file)
  try(system(destDirRSEMCount.create))
  destDirRSEMCount
}
create.Cufflinks.folder <- function (dest.dir, text.add,log.file=NULL){
  if(is.null(dest.dir) || is.null(text.add)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  destDirCufflinksCount <-  paste(dest.dir, text.add, ".Cufflinks.Counting/", sep="")
  destDirCufflinksCount.create <- paste("mkdir ", destDirCufflinksCount, sep="")
  #
  add.to.logs(destDirCufflinksCount.create,log.file)
  try(system(destDirCufflinksCount.create))
  destDirCufflinksCount
}
create.Collect.folder <- function (dest.dir, text.add,log.file=NULL){
  if(is.null(dest.dir) || is.null(text.add)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  collectDir <- paste(dest.dir,text.add,".Collect/",sep="")
  collectDir.create <- paste("mkdir",collectDir)
  
  add.to.logs(collectDir.create,log.file)
  try(system(collectDir.create))
  collectDir
}
create.oldRun.Collect.folder <- function(previous.dir,text.add,log.file=NULL){
  if(is.null(previous.dir) || is.null(text.add)) stop("Arguments should not be NULL")
  previous.dir <- trailDirCheck(previous.dir)
  collectDir <- paste(previous.dir,text.add,".Collect/",sep="")
  collectDir.create <- paste("mkdir",collectDir)
  add.to.logs(collectDir.create,log.file)
  try(system(collectDir.create))
  collectDir
}
copy.attribute.file.to.dest <- function(attrPath,dest.dir,log.file=NULL){
  if(is.null(dest.dir) || is.null(attrPath)) stop("Arguments should not be NULL")
  add.to.logs("################## Copying attribute file for this run into destination folder ##################",log.file)
  dest.dir <- trailDirCheck(dest.dir)
  copyAttributeFile <- paste("cp",attrPath,dest.dir)
  add.to.logs(copyAttributeFile,log.file)
  try(system(copyAttributeFile))
  # Returns command for unit tests
  copyAttributeFile
}
create.run.logs <- function(destDirLog,text.add) {
  ###############
  # Log-tail command addition
  ###############
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
  for (i in 1: length(fqfiles.base)) {
    # This treats all files like .fq file
    left.fq <- paste(fqfiles.base[i],".1.fq", sep="")
    right.fq <- paste(fqfiles.base[i],".2.fq", sep="")
    fqfiles.table <- rbind(fqfiles.table, c(left.fq, right.fq)) }
  fqfiles.table
}
##########################################################################
########################## NO DATA PREP ##################################
##########################################################################
rename.preprocessed.files <- function(dest.dir){
  # 
  # Rename .fastq files to .fq files
  #
  if (is.null(dest.dir)) stop("Arguments should not be NULL")
  setwd(dest.dir)
  active.files <- dir(dest.dir)
  rename.command <- paste("date")
  for (num in 1:length(active.files)){
    if (grepl(".fastq$",active.files[num])){
      to.rename.base <- substr(active.files[num],1,nchar(active.files[num])-5)
      move.file <- paste(to.rename.base,"fq",sep="")
      rename.command <- paste(rename.command,";","mv",active.files[num],move.file)
    } 
  }
  system(rename.command)
  # Returns command for unit tests
  rename.command
}
# Copying the ready-to-process fastq files to the destination directory:
copy.preprocessed.files <- function(readyData.dir,dest.dir,log.file=NULL) {
  if (is.null(readyData.dir) || is.null(dest.dir)) stop("Arguments should not be NULL")
  readyData.dir <- trailDirCheck(readyData.dir)
  dest.dir <- trailDirCheck(dest.dir)
  copyFq <- paste(readyData.dir,"*.fq",sep="")
  readyDataCopy <- paste("cp",copyFq,dest.dir)
  active.files <- dir(dest.dir)
  fastq.files <- FALSE
  for (num in 1:length(active.files)){
    if (grepl(".fastq$",active.files[num])){
      copyFastq <- paste(readyData.dir,"*.fastq",sep="")
      readyDataCopy <- paste(readyDataCopy,";",copyFastq)
      fastq.files <- TRUE
    }
  }
  add.to.logs("################## Copying preprocessed .fq files into destination directory  ##################",log.file)
  add.to.logs(readyDataCopy,log.file)
  try(system(readyDataCopy))
  if (fastq.files){
    rename.preprocessed.files(dest.dir)
  }
  # Returns command for unit tests
  readyDataCopy
}
convert.file.list.to.table <- function(paired.end,dest.dir) {
  if (is.null(paired.end) || is.null(dest.dir)) stop("Arguments should not be NULL")
  # Paired End
  if (paired.end) {
    fqfiles <- dir(dest.dir) 
    fqfiles <- fqfiles[grep(".fq$", fqfiles)]
    fqfiles.base <- substr(fqfiles, 1,nchar(fqfiles) - 5)
    fqfiles.base <- unique(fqfiles.base)
    # Table of pairs of FASTQ file names (1 pair per row):
    # We assemble left and right file names explicitly here because 
    # We don't want to rely on file sequence in the directory and get wrong results if there are occasional unpaired files etc. 
    # Assemble function for paired-end libraries:
    fqfiles.table <- fqfiles.table.pe.assemble(fqfiles.base)
  }
  # Single End
  if (!(paired.end)) {
    fqfiles <- dir(dest.dir)
    fqfiles <- fqfiles[grep(".fq$", fqfiles)]
    fqfiles.table <- cbind(NULL, fqfiles)
  }
  fqfiles.table
}
##########################################################################
########################## DATA PREP #####################################
##########################################################################
# If the data needs to be processed all the way
# from the compressed archives in the raw directory
find.files.for.dataprep <- function(raw.dir,unzipped,log.file=NULL){
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
      add.to.logs("################## Finding zipped files to prepare ##################",log.file)
      gzfiles <- gzfiles[grep(".gz", gzfiles)]
      if(!presplit) fqfiles.base <- substr(gzfiles, 1,nchar(gzfiles) - 3)
      if (presplit) fqfiles.base <- substr(gzfiles, 1,nchar(gzfiles) - 3)
      fqfiles.table <- fqfiles.table.pe.assemble(fqfiles.base)
    }
    if (unzipped) {
      add.to.logs("################## Finding unzipped files to prepare ##################",log.file)
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
      add.to.logs("################## Finding zipped files to prepare ##################",log.file)
      gzfiles <- dir(raw.dir) 
      gzfiles <- gzfiles[grep(".gz", gzfiles)] # just in case raw.dir has somethig else besides .gz files
      fqfiles <- substr(gzfiles, 1, nchar(gzfiles)-3)
      fqfiles.table <- cbind(NULL, fqfiles)
    }
    if (unzipped){
      add.to.logs("################## Finding unzipped files to prepare ##################",log.file)
      fqfiles <- dir(raw.dir) 
      fqfiles <- fqfiles[grep(".fq|.fastq", fqfiles)] # just in case raw.dir has somethig else besides .fq or .fastq files
      fqfiles <- substr(fqfiles, 1, nchar(fqfiles)-0)
      fqfiles.table <- cbind(NULL, fqfiles)
    }
    fqfiles.table
  }
}
##########################################################################
########################## RANGES OF DATA ################################
##########################################################################
prepare.chunk.function <- function(fqfiles.table,nStreams){
  if (is.null(fqfiles.table) || is.null(nStreams)) stop("Arguments should not be NULL")
  nStreams <- check.nStreams(fqfiles.table,nStreams)
  chunk <- function(x, n) split(x, sort(rank(x) %% n)) # commonly known solution to divide data equally
  rangelist <- chunk(1:nrow(fqfiles.table), nStreams)
  # ranges are in rangelist[1:nSreams]
  for (ii in 1:nStreams) assign(paste("range",ii,sep=""),rangelist[[ii]])
  rangelist
}

check.nStreams <- function(fqfiles.table,nStreams){
  if (nStreams > nrow(fqfiles.table)) nStreams <- nrow(fqfiles.table) # quiet and nice solution but there will de descrepancy between nStreams in the attribute file and the actual nStreams
  nStreams
}
##########################################################################
###################### EXPRESSION QUANTIFICATION #########################
##########################################################################
# Some alignments have special case commands which need to be handled
# This will allow for that to be done in the top script
#
# This special case is to facilitate CUSHAW-GPU which processes data
# A bit differently then all the other alignment protocols
# Doesn't run the Alignment if no expresion calculation is desired
start.alignment.process <- function(base.dir,rangelist,nStreams,log.file=NULL){
  add.to.logs(rangelist,log.file)
  setwd(base.dir)
  GPUspecialCase <- FALSE
  this.resName <- NULL
  add.to.logs("################## Creating alignment script ##################",log.file)
  source("GLSeq.Alignment.R")
  # RETURNS THE COMM STACK
  comm.stack.pool
}
start.counting.process <- function(countable.sams.dir,dest.dir,base.dir,log.file=NULL){
  add.to.logs(paste("Starting counting:", fqfiles.table),log.file)
  #
  comm.stack.pool <- "date"
  #
  countable.sams <- dir(countable.sams.dir)[grep("countable.sam", dir(countable.sams.dir))]
  #
  for (i in 1:length(countable.sams)){
    countable.sam <- countable.sams[i]
    copy.comm <- paste("cp",paste(countable.sams.dir,countable.sam,sep="/"),dest.dir)
    #
    add.to.logs("################## Copying premade SAM files into directory ##################",log.file)
    add.to.logs(copy.comm,log.file)
    #
    system(copy.comm)
    setwd(base.dir)
    source("GLSeq.Counting.R")
  }
  comm.stack.pool <- paste(comm.stack.pool,"&&",count.comm,"&")
  
  # RETURNS THE COMM STACK
  comm.stack.pool
}
#
RSEM.finish <- function(comm.stack.pool,destDirRSEMCount,dest.dir) {
  if (is.null(comm.stack.pool) || is.null(destDirRSEMCount) || is.null(dest.dir)) stop("Arguments should not be NULL")
    destDirRSEMCount <- trailDirCheck(destDirRSEMCount)
    dest.dir <- trailDirCheck(dest.dir)
    comm.stack.pool <- paste(comm.stack.pool,"wait")
    comm.stack.pool <- paste(comm.stack.pool,"&&","cd",dest.dir,"&&","mv",paste("*.RSEM.counts.*"),destDirRSEMCount)
    comm.stack.pool <- paste(comm.stack.pool,"&&","mv","*.index.*",destDirRSEMCount)   
    # RETURNS THE COMM STACK
    comm.stack.pool
}
##########################################################################
####################### SAVE RUN VARIABLES ###############################
##########################################################################
save.run.data <- function(base.dir,text.add) {
  setwd(base.dir)
  currentRun.dataFile <- paste("GLSeq.vars.", text.add, ".rda", sep="")
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
#
##########################################################################
####################### PREPARE DATA #####################################
##########################################################################
run.data.prep <- function(destDirLog,text.add,attrPath,dest.dir){
  if (is.null(destDirLog) || is.null(text.add) || is.null(attrPath) || is.null(dest.dir)) stop("Arguments should not be NULL")
  destDirLog <- trailDirCheck(destDirLog)
  dest.dir <- trailDirCheck(dest.dir) 
  Sys.sleep(10) 
  dataPrepLog <- paste(destDirLog, text.add, ".DataPrepLog.txt", sep="")
  dataPrepErr <- paste(destDirLog, text.add, ".DataPrepErr.txt", sep="")
  dataPrep <- paste("Rscript GLSeq.dataprep.R ", text.add," ",dest.dir," ",attrPath," 1>> ", dataPrepLog, " 2>> ", dataPrepErr, sep="") 
  print("Starting Data Preparation")
  system(dataPrep)
  #
  #################
  # Watch for readiness of the data files 
  # before starting the expression calculations
  ################
  DataIsWaiting <- TRUE
  if (dataPrepare == "dataprep"){
    setwd(dest.dir)
    DataIsWaiting <- FALSE 
  }
  dataReady.ind <- paste(text.add, ".DataReady", sep="") 
  # Waiting
  while(!(DataIsWaiting)) {
    Sys.sleep(21)
    DataIsWaiting <- dataReady.ind %in% dir(dest.dir)
  }
  if (DataIsWaiting) print("Data Preparation Complete")
}
##########################################################################
#################### EXPRESSION CALCULATION ##############################
##########################################################################
cushaw.gpu.run <- function (Cushawgpu.special.case,log.file=NULL) {
  if (is.null(Cushawgpu.special.case)) stop("Arguments should not be NULL")
  add.to.logs("################## Alignment with CUSHAW-GPU starting ##################",log.file)
  add.to.logs(Cushawgpu.special.case,log.file)
  cushaw.start.time <- proc.time()
  system(Cushawgpu.special.case)
  #
  add.to.logs(paste("Cushaw alignment took", proc.time()[3] - cushaw.start.time[3],"seconds to complete"),log.file)
  #
  warning("The alignment step has now completed.")
}
execute.comm.stack <- function (comm.stack.pool, log.file=NULL) {
  if (is.null(comm.stack.pool)) stop("Arguments should not be NULL")
  add.to.logs("################## Executing communication script ##################",log.file)
  add.to.logs(comm.stack.pool,log.file)
  stack.start.time <- proc.time()
  t <- try(system(comm.stack.pool))
  #
  add.to.logs(paste("Alignment and/or counting stack took",proc.time()[3] - stack.start.time[3],"seconds to complete"),log.file)
  #
  if (t == 0){
    add.to.logs("Communication script executed successfully.",log.file)
    print("Successfully executed command")
  } else{
    add.to.logs("Communication script exited with a failure message.",log.file)
  }
}
#
###########################################################################
###################### COLLECT ############################################
###########################################################################
#
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
  add.to.logs("################## Starting Results Collection ##################",log.file)
  print("STARTED COLLECTION")
  source("GLSeq.ResultsCollect.R")
}