#########################################################
# Great Lakes Seq package for low-level processing of RNA-Seq data
# Trimmomatic-BWA-HTSeq quantification of the expression values
# (paired end libraries); based on "Dana's script" + Irene input + Oleg's observations on the data
# (results of the latter: default headcrop is increased to 12 and Illumina artifacts removal is introduced)
# September 2013
#########################################################
#
#
source("GLSeq.Util.R")
# if (!(paired.end)) stop('BWA pipeline is currently supported for paired-end libraries only \n')
setwd(dest.dir)
# copy genome indices to the destimation dir: 
ref.dir <- paste(base.dir, rGenome, sep="")
indCopy <- paste("cd ", ref.dir, " && cp * ", dest.dir, sep="")
system(indCopy)
#
comm.stack.pool <- NULL # 
#
for (zz in 1:nStreams) {
# assembly and runing the system command, one library at a time:
for (i in rangelist[[zz]]) {
###################
# Alignment with SAM output
###################
# names of current fastq files:
fq.left <- fqfiles.table[i,1]
if (paired.end) fq.right <- fqfiles.table[i,2]
this.library <- substr(fqfiles.table[i,1], 1,libNchar)
this.resName <- paste(this.library, text.add, sep=".")
unsorted.sam <- paste(this.resName, "unsorted.sam", sep=".")
#
###################
# CUSHAW alignment:
################### 
	if (paired.end) {
if (!(GPU.accel)) sam.create <- paste(CUSHAW.path, "-r", refFASTAname, "-q", fq.left, fq.right, "-o", unsorted.sam, "-rgpl", seqPlatform, "-t", nCores) 
if (GPU.accel) sam.create <- paste(CUSHAW.GPU.path, "-r", refFASTAname, "-q", fq.left, fq.right, "-o", unsorted.sam, "-rgpl", seqPlatform, "-t", nCores) }
	if (!(paired.end)) {
if (!(GPU.accel)) sam.create <- paste(CUSHAW.path, "-r", refFASTAname, "-s", fq.left, "-o", unsorted.sam, "-rgpl", seqPlatform, "-t", nCores) 
if (GPU.accel) sam.create <- paste(CUSHAW.GPU.path, "-r", refFASTAname, "-s", fq.left, "-o", unsorted.sam, "-rgpl", seqPlatform, "-t", nCores) }
#
###################
# SAM => BAM file with index 
###################
sorted.arg <- paste(this.resName, "sorted", sep=".")
ref.index <- paste(refFASTAname, "fai", sep=".") 
sorted.bam <- paste(this.resName, "sorted.bam", sep=".") 
bam.create <- paste("samtools view -uS -t ", ref.index, unsorted.sam, " | samtools sort - ", sorted.arg) #
bam.index <- paste("samtools index", sorted.bam) # 
#
###################
# Converting the final bam (coordinate-sorted)
# to name-sorted  sam:  
###################
#
countable.sam <- paste(this.resName, "countable.sam", sep=".")
paired.arg <- paste(this.resName, "paired", sep=".")
paired.bam <- paste(this.resName, "paired.bam", sep=".")
# let's avoid pipe here to save RAM
# Directly write the samtools output to a new countable file after sorting
# -n in the samtools sort indicates that the document will be sorted by read name rather than chromosome
# -h in the samtools view indicates a header will be printed
countable.comm <- paste("samtools sort -n", sorted.bam, paired.arg, "&&", "samtools view -h", paired.bam, ">", countable.sam) # System command #8
###################
# Current command:
###################
# current command: 
if (paired.end) comm.i <- paste(sam.create, "&&", bam.create, "&&", bam.index)
if (!(paired.end)) comm.i <- paste(sam.create, "&&", bam.create, "&&", bam.index)
# for the very first assembly in the stack: 
if (i == rangelist[[zz]][1])  comm.stack.pool <- paste(comm.stack.pool, " date && ", comm.i)
# for subsequent assemblies of every stack: 
if (i != rangelist[[zz]][1])  comm.stack.pool <- paste(comm.stack.pool, " && date && ", comm.i)
# system(comm.i)
} # for i 
if (zz ==1) fileCompletenessID <- paste(text.add, ".completeExpression", sep="")
comm.stack.pool <- paste(comm.stack.pool,  " && echo  >  ", fileCompletenessID, ".", zz, " & ", sep="")
} # for zz
#
# collection of the results: 
collLog <- paste(destDirLog, text.add, ".ResultsCollectLog.txt", sep="")
collerr <- paste(destDirLog, text.add, ".ResultsCollectErrors.txt", sep="")
collResults <- paste("cd ", base.dir, " && ", "Rscript GLSeqResultsCollect.R ", text.add, base.dir, dest.dir, " 0 1>> ", collLog, " 2>> ", collerr, " &", sep="") 
if (resCollect == "nocollect") collResults <- "\n"
# 
# pool of all the system commands (versions for +/- compute expression and +/- collect results): 
#
comm.stack.pool <- paste(comm.stack.pool, collResults, "\n", sep=" ") 
if (exprRun == "noexprcalc") comm.stack.pool <- paste(collResults, "\n", sep=" ") 
#
###########################################################################




