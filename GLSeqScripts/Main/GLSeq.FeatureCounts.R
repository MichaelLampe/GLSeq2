args <- commandArgs(trailingOnly = TRUE)
countable.sam <- as.character(args[1])
destDirFeatureCounts <- as.character(args[2])
refGFFname <- as.character(args[3])
dest.dir <- as.character(args[4])
this.resName <- as.character(args[5])
paired.end <- as.logical(args[6])
idAttr <- as.character(args[7])
FeatureCountsSpecialOptions <- as.character(args[8])

library(Rsubread)
#############
# Makes sure in the correct dir
#############
setwd(dest.dir)
#
#############
# Reads output into a file
#############
counts.summary <- paste(this.resName,".FeatureCounts.summary.txt",sep="")
counts.data <- paste(this.resName,".FeatureCounts.counts.csv",sep="")
counts.stats <- paste(this.resName,".FeatureCounts.stats.csv",sep="")
counts.annotation <- paste(this.resName,".FeatureCounts.annotations.csv",sep="")
# If we've been given a junk countable.sam name, ignore it.
if (!is.na(countable.sam)) {
  sink(file = counts.summary, append = TRUE, type = c("output"), split = FALSE)
  # Default advanced feature options
  requireBothEndsMapped <- FALSE
  countChimericFragments <- FALSE
  minFragLength <- 50
  maxFragLength <- 600
  minMQS <- 0
  autosort <- FALSE
  largestOverlap <- FALSE
  minOverlap <- 1
  countPrimaryAlignmentsOnly <- FALSE
  readExtension5 <- 0
  readExtension3 <- 0
  # by splitting on -, each will correspond to a command
  command <- unlist(strsplit(FeatureCountsSpecialOptions, "-"))
  # Iterate through the commands and parse them for values
  # If command is 0 (Aka no command, it will try to make a list 1,0, so we need to exclude that)
  if (length(command) != 0){
    for (commandIndex in 1:length(command)) {
      # Trim whitespace on the front and end
      if (!is.na(command[commandIndex])){
        value <- gsub("^\\s+|\\s+$", "",command[commandIndex])
        if (command[commandIndex] == "B") {
          requireBothEndsMapped <- TRUE
        }
        else if (command[commandIndex] == "C") {
          countChimericFragments <- TRUE
        }
        else if (grepl("d", command[commandIndex])){
          # Chop off the command and just take the value
          # Then convert the string to an int
          value <- unlist(strsplit(command[commandIndex],"d "))[2]
          minFragLength <- strtoi(value)
        }
        else if (grepl("D", command[commandIndex])){
          # Chop off the command and just take the value
          # Then convert the string to an int
          value <- unlist(strsplit(command[commandIndex],"D "))[2]
          maxFragLength <- strtoi(value)
        }
        else if (grepl("Q", command[commandIndex])){
          # Chop off the command and just take the value
          # Then convert the string to an int
          value <- unlist(strsplit(command[commandIndex],"Q "))[2]
          minMQS <- strtoi(value)
        }
        else if ("donotsort" == command[commandIndex]) {
          autosort <- TRUE
        }
        else if ("largestOverlap" == command[commandIndex]) {
          largestOverlap <- TRUE
        }
        else if (grepl("minOverlap", command[commandIndex])) {
          value <- unlist(strsplit(command[commandIndex],"minOverlap "))[2]
          minOverlap <- strtoi(value)
        }
        else if ("primary" == command[commandIndex]){
          countPrimaryAlignmentsOnly <- TRUE
        }
        else if (grepl("readExtension5", command[commandIndex])) {
          value <- unlist(strsplit(command[commandIndex],"readExtension5 "))[2]
          readExtension5 <- strtoi(value)
        }
        else if (grepl("readExtension3", command[commandIndex])) {
          value <- unlist(strsplit(command[commandIndex],"readExtension3 "))[2]
          readExtension3 <- strtoi(value)
        }
      }
    }
  }

  fc <- featureCounts(
    files=countable.sam,
    annot.ext=paste(dest.dir,refGFFname,sep=""),
    isGTFAnnotationFile=TRUE,
    GTF.featureType="exon",
    GTF.attrType="gene_id",
    isPairedEnd=paired.end,
    nthreads=1,
    requireBothEndsMapped=requireBothEndsMapped, # Start of advanced features
    countChimericFragments=countChimericFragments,
    minFragLength=minFragLength,
    maxFragLength=maxFragLength,
    minMQS=minMQS,
    autosort=autosort,
    largestOverlap=largestOverlap,
    minOverlap=minOverlap,
    countPrimaryAlignmentsOnly=countPrimaryAlignmentsOnly,
    readExtension5=readExtension5,
    readExtension3=readExtension3
  )
  sink()
  write.csv(fc[1],counts.data)
  write.csv(fc[2],counts.annotation)
  write.csv(fc[4],counts.stats)
}