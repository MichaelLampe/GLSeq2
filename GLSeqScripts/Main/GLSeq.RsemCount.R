source ("GLSeq.Util.R")
#############
# Null Fixes
#############
#
if (is.null(this.resName)) this.resName <- paste(dest.dir,text.add,sep="")
#
#
# Prints out a check to the log file if the counting file is correct.
#
header.message <- paste("echo \"Checking SAM files for RSEM compatability. The file is: ",countable.sam,"\"",sep="")
destDirRSEMCount <- trailDirCheck(destDirRSEMCount)
rsem.log.file <- paste(destDirRSEMCount,fqfiles.table[i,1],"_File_Validation",sep="")
header.file <- paste(header.message,">> ",rsem.log.file)
check.file <- paste("rsem-sam-validator",countable.sam,">>",rsem.log.file)
#
refName <- paste(this.resName,"RSEM","index",sep=".")
prepareReference <- paste("rsem-prepare-reference",paste(dest.dir,refFASTAname,sep=""),refName)
#
countFile <- paste(this.resName,"RSEM","counts",sep=".")
if (Condor){
  # While we set samtools to having a lot of memory to use, it will rarely use that much.
  # Because of how flexible condor is, we can set it to only 6GB and then it will scale up if needed
  rsemOptions <- paste("--calc-pme","--calc-ci","-p",8,"--samtools-sort-mem 32G")
} else{
  rsemOptions <- paste("--calc-pme","--calc-ci","-p",nCores)
}
if (paired.end) rsemOptions <- paste(rsemOptions,"--sam","--paired-end")
if (!paired.end) rsemOptions <- paste(rsemOptions,"--sam")
# Organize the files into the RSEM count folder.
organize.files <- paste("mv",paste(this.resName,"*.RSEM.counts.*",sep=""),destDirRSEMCount)
organize.files <- paste(organize.files,"&&","mv",paste(this.resName,"*.index.*",sep=""),destDirRSEMCount)

if (!is.na(countable.sam)){
  calculateExpression <- paste("rsem-calculate-expression",rsemOptions,countable.sam,refName,countFile)
  if (count.comm != "") count.comm <- paste(count.comm,"&&",check.file,"&&", prepareReference, "&&",calculateExpression,"&&",organize.files)
  if (count.comm == "") count.comm <- paste(header.file,"&&",check.file,"&&",prepareReference, "&&", calculateExpression,"&&",organize.files)
}