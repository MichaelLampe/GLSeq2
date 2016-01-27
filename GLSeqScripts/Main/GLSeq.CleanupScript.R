# Create cleanup directory
cleanup.directory <- paste(trailDirCheck(storage.destination),text.add,sep="")
create.cleanup.directory <- paste("mkdir",cleanup.directory)
executeOrPrint(create.cleanup.directory, Condor)

# Move Countable SAM files and results into that file
move.files <- "mv *.Counting & mv *.Collect & mv *.countable.sam"
executeOrPrint(move.files, Condor)
