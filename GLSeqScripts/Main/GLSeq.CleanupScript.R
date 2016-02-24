source("GLSeq.Util.R")
# Create cleanup directory
cleanup.directory <- paste(trailDirCheck(storage.destination),text.add,sep="")
create.cleanup.directory <- paste("mkdir",cleanup.directory)
printOrExecute(create.cleanup.directory, Condor)

# Move Countable SAM files and results into that file
move.files <- "cp *.Counting & cp *.Collect & cp *.countable.sam"
printOrExecute(move.files, Condor)
