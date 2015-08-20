args <- commandArgs(trailingOnly = TRUE)
load.file <- as.character(args[1])
load(load.file)

setwd(destDirHTSeqCount)
#
###########################
# Generating normalized counts (FPKM) and saving the file to disk:
# function to pull GeneID from the first record of a particular element of the long character string ("mess") in the column 9 of the gtf file:
###########################
genePull <- function(gtf.gene.record) {
  geneID.record <- strsplit(gtf.gene.record, split=";")[[1]][1]
  geneID <- substr(geneID.record, 9, nchar(geneID.record))
  geneID
}
#
##########################
# Extracting an arbitrary IDs from the annotation column of a gtf/gff file:
##########################
idPull <- function(annotColumn, ID2look) {
  annotColumn.split <- strsplit(annotColumn, split=";")
  ID2look.positions <- lapply(annotColumn.split, grep, pattern=ID2look)
  IDvec <- rep(NA, length(annotColumn.split))
  for (vv in 1:length(IDvec)) {
    if (length(ID2look.positions[[vv]]) > 0) {
      IDcopyNum <- length(ID2look.positions[[vv]]) # the LAST copy of the same ID is actually used in the count computation!
      IDvec[vv] <- annotColumn.split[[vv]][ID2look.positions[[vv]]][IDcopyNum]
      IDvec[vv] <- substr(IDvec[vv], nchar(ID2look)+2, nchar(IDvec[vv]))
    }
  } # for vv
  # The length of the original annotColumn is preserved (records not containing the ID2look now contain NAs)
  IDvec
}
#
#
cfiles <- dir(pattern="counts")
countDirName <- paste(text.add, "counts", sep=".")
cfiles <- cfiles[cfiles != countDirName]
count.mtrx <- NULL
#
for (i in 1:length(cfiles)) {
  data.i <- read.table(cfiles[i], row.names=1)
  colnames(data.i) <- substr(cfiles[i],1,libNchar)
  if (is.null(count.mtrx)) count.mtrx <- data.i
  if (!(is.null(count.mtrx)) & i != 1) count.mtrx <- cbind(count.mtrx, data.i)
}
########################################
# Generating file with report on all abnormal events during alignment:
########################################
#
# Create the Exception Report CSV
# drop=F is required to keep the data as a data frame
# which allows for names to be retained
#
data.exceptions <- data.i[(nrow(data.i)-4):nrow(data.i),,drop=F]
exrep.fName <-  paste(collectDir, text.add, ".HTSeq.exceptionReport.csv", sep="")
write.csv(data.exceptions, file=exrep.fName)
########################################
# Generating file with report on all other events during alignment:
########################################
# Create the Counts Report CSV
# drop=F is required to keep the data as a data frame
# which allows for names to be retained
#
count.mtrx <- count.mtrx[1:(nrow(data.i)-5),,drop=F]
HTSeq.counts.fName <-  paste(collectDir, text.add, ".HTSeq.counts.csv", sep="")
write.csv(count.mtrx, file=HTSeq.counts.fName)
#
# Load in the GTF file
#
reference.gtf <- paste(dest.dir,refGFFname,sep="")
gtf <-  read.table(reference.gtf, sep="\t", header=FALSE, as.is=TRUE)
mess <- gtf[,gtfFeatureColumn] # 9
geneIDs.full <- idPull(mess, idAttr)
lengthData <- cbind(gtf[,3:5], geneIDs.full) # columns: 1) feature, 2) start, 3) end, 4) ID
lengthData <- lengthData[lengthData[,1] == "exon",]
#
# we may see situations when "exon" record does not have  desireable ID record in the gtfFeatureColumn:
#
lengthData <- lengthData[!(is.na(lengthData[,4])),]
lengthData <- cbind(lengthData, lengthData[,3] - lengthData[,2])
colnames(lengthData)[5] <- "length"
#
# Considering multiple exons per gene:
#
duplicated.IDs.summary <- table(lengthData[,4])[table(lengthData[,4]) > 1]
lengthData.singleExons <- lengthData[!(lengthData[,4] %in% names(duplicated.IDs.summary)),]
lengthData.multipleExons <- lengthData[lengthData[,4] %in% names(duplicated.IDs.summary),]
lengthData <- lengthData.singleExons # will be extended to summarized exon length below
for (multiexon in 1:length(names(duplicated.IDs.summary))) {
  data.multiexon <- lengthData.multipleExons[lengthData.multipleExons[,4] == names(duplicated.IDs.summary)[multiexon],]
  #
  # Recording the sum of all exon's lengths for a multi-exon gene in the first line of the length data for this gene:
  #
  data.multiexon[1,5] <- sum(data.multiexon[,5])
  lengthData <- rbind(lengthData, data.multiexon[1,])
} # for multiexon
#
# End of multiple exons per gene summarization
#
lengthData <- lengthData[lengthData[,4] %in% rownames(count.mtrx),] # avoiding trouble with occasional NAs in the column 4 of the lengthData
rownames(lengthData) <- lengthData[,4]
#
# Sorting / restricting the count matrix based on the non-redundant ID vector:
#
lengthData <- lengthData[rownames(count.mtrx),]
#
for (normCol in 1:ncol(count.mtrx)) {
  mil.mappedR <- sum(count.mtrx[,normCol]) / 1000000
  kbase.transcr <- lengthData[,5] / 1000
  count.mtrx[,normCol] <- count.mtrx[,normCol] / (kbase.transcr * mil.mappedR)
} # for normCol
#
# Create the normalized RPKM CSV
#
HTSeq.RPKM.fName <-  paste(collectDir, text.add, ".HTSeq.RPKM.csv", sep="")
write.csv(count.mtrx, file=HTSeq.RPKM.fName)
#
# Visualization files (the name is common with RSEM to unify the code for strandExtract version:
#
bam.files.thisRun.genome <- dir(pattern="*.bam")
bai.ind <- grep("bam.bai", bam.files.thisRun.genome)
bam.files.thisRun.genome <- bam.files.thisRun.genome[-bai.ind]
bai.files.thisRun.genome <- bam.files.thisRun.genome[bai.ind]
#
#
if(!(strandExtract)) {
  vizfiles.base <- substr(bam.files.thisRun.genome, 1, nchar(bam.files.thisRun.genome)-4)
  wigfiles <- paste(vizfiles.base, "wig", sep=".")
  bamMove <- paste("mv *.bam*", destDirBam)
  for (ii in 1:length(bam.files.thisRun.genome)) {
    # wigGen.ii <- paste("rsem-bam2wig",  bamfiles[ii], wigfiles[ii], vizfiles.base[ii])
    wigGen.ii <- paste(bam2wigPath, bam.files.thisRun.genome[ii], wigfiles[ii])
    if (ii == 1) wigGen <- wigGen.ii
    if (ii != 1) wigGen <- paste(wigGen, "&", wigGen.ii)
  } # for bamfiles
  try(system(bamMove))
  Sys.sleep(5)
  setwd(destDirBam)
  system(paste(wigGen, "&"))
  Sys.sleep(2)
} # If not strand extract
#
######################################################
################## PREPARE GRAPHS ####################
######################################################
#
setwd(collectDir)
HTSeq.counts <- read.csv(HTSeq.counts.fName, header = TRUE, row.names=1)
HTSeq.FPKM <- read.csv(HTSeq.RPKM.fName, header = TRUE, row.names=1)
#
HTSeq.counts.graphTitle <- "HTSeq log10(Counts+1)"
HTSeq.FPKM.graphTitle <- "HTSeq log10(FPKM+0.01)"
graph.title <- paste(text.add,".HTSeq.BoxPlot.png",sep="")
png(graph.title,width=1400,height=700)
par(mfrow=c(1,2))
boxplot(log10(HTSeq.counts + 1), main=HTSeq.counts.graphTitle)
boxplot(log10(HTSeq.FPKM + 0.01), main=HTSeq.FPKM.graphTitle)
dev.off()