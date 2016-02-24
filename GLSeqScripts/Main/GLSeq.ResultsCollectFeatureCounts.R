args <- commandArgs(trailingOnly = TRUE)
load.file <- as.character(args[1])
load(load.file)
setwd(destDirFeatureCountsCount)
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
      IDcopyNum <- length(ID2look.positions[[vv]]) # the LAST copy of the same ID is actullay used in the count computation!
      IDvec[vv] <- annotColumn.split[[vv]][ID2look.positions[[vv]]][IDcopyNum]
      IDvec[vv] <- substr(IDvec[vv], nchar(ID2look)+2, nchar(IDvec[vv]))
    }
  } # for vv
  # The length of the original annotColumn is preserved (records not containing the ID2look now contain NAs)
  IDvec
}
#
cfiles <- dir(pattern=".FeatureCounts.counts.csv")
countDirName <- paste(text.add, "counts", sep=".")
cfiles <- cfiles[cfiles != countDirName]
count.mtrx <- NULL
#
#
# Need to fix the formatting of the file because FeatureCounts has a bunch of useless/artistic outputs that clutter the file.
#
# correctedNames <- c()
# for (i in 1:length(cfiles)){
#   setwd(destDirFeatureCountsCount)
#   correctedName = paste(cfiles[i],".corrected",sep="")
#   correctedNames <- c(correctedNames,correctedName)
#   openFile <- file(cfiles[i], open="r")
#   document <- readLines(openFile)
#   for (i in 1:length(document))
#     if (grepl("$counts",document[i],fixed=TRUE)){
#       i <- i + 2
#       setwd(collectDir)
#       while (!(grepl("$annotation",document[i],fixed=TRUE))){
#         write(document[i],file=correctedName,append=TRUE)
#         i <- i + 1
#       }
#       currentLine <- readLines(openFile)
#     }
#   close(openFile)
# }
# cfiles <- correctedNames
#
setwd(collectDir)
for (i in 1:length(cfiles)) {
  data.i <- read.csv(paste(destDirFeatureCountsCount,cfiles[i],sep=""),row.names=1)
  #data.i <- read.table(cfiles[i], row.names=1)
  colnames(data.i) <- cfiles[i]
  if (is.null(count.mtrx)) count.mtrx <- data.i
  if (!(is.null(count.mtrx)) & i != 1) count.mtrx <- cbind(count.mtrx, data.i)
}
#
featureCounts.counts.fName <-  paste(collectDir, text.add, ".FeatureCounts.counts.csv", sep="")
write.csv(count.mtrx, file=featureCounts.counts.fName)
#
#
# Load GTF file
reference.gtf <- paste(dest.dir,refGFFname,sep="")
gtf <-  read.table(reference.gtf, sep="\t", header=FALSE, as.is=TRUE)
mess <- gtf[,9] # 9
geneIDs.full <- idPull(mess, idAttr)
lengthData <- cbind(gtf[,3:5], geneIDs.full) # columns: 1) feature, 2) start, 3) end, 4) ID
lengthData <- lengthData[lengthData[,1] == "exon",]
# we may see situations when "exon" record does not have  desireable ID record in the gtfFeatureColumn:
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
  # recording the sum of all exon's lengths for a multi-exon gene in the first line of the length data for this gene:
  data.multiexon[1,5] <- sum(data.multiexon[,5])
  lengthData <- rbind(lengthData, data.multiexon[1,])
} # for multiexon


# End of multiple exons per gene summarization
#
lengthData <- lengthData[lengthData[,4] %in% rownames(count.mtrx),] # avoiding trouble with occasional NAs in the column 4 of the lengthData
rownames(lengthData) <- lengthData[,4]

# Sorting / restricting the count matrix based on the non-redundant ID vector:
lengthData <- lengthData[rownames(count.mtrx),]
for (normCol in 1:ncol(count.mtrx)) {
  mil.mappedR <- sum(count.mtrx[,normCol]) / 1000000
  kbase.transcr <- lengthData[,5] / 1000
  count.mtrx[,normCol] <- count.mtrx[,normCol] / (kbase.transcr * mil.mappedR)
} # for normCol
featureCounts.RPKM.fName <-  paste(collectDir, text.add, ".FeatureCounts.RPKM.csv", sep="")
write.csv(count.mtrx, file=featureCounts.RPKM.fName)
#
# visualization files (the name is common with RSEM to unify the code for strandExtract version:
bam.files.thisRun.genome <- dir(pattern="*.bam")
bai.ind <- grep("bam.bai", bam.files.thisRun.genome)
bam.files.thisRun.genome <- bam.files.thisRun.genome[-bai.ind]
bai.files.thisRun.genome <- bam.files.thisRun.genome[bai.ind]
#
#





# Going to temp disable this for now.

# Permission denied error (That's on Oleg's end I think)

#if(FALSE) {
 # vizfiles.base <- substr(bam.files.thisRun.genome, 1, nchar(bam.files.thisRun.genome)-4)
 # wigfiles <- paste(vizfiles.base, "wig", sep=".")
 # bamMove <- paste("mv *.bam*", destDirBam)
 # for (ii in 1:length(bam.files.thisRun.genome)) {
    # wigGen.ii <- paste("rsem-bam2wig",  bamfiles[ii], wigfiles[ii], vizfiles.base[ii])
 #   wigGen.ii <- paste(bam2wigPath, bam.files.thisRun.genome[ii], wigfiles[ii])
 #   if (ii == 1) wigGen <- wigGen.ii
 #   if (ii != 1) wigGen <- paste(wigGen, "&", wigGen.ii)
 # } # for bamfiles
 # try(system(bamMove))
 # Sys.sleep(5)
 # setwd(destDirBam)
 # system(paste(wigGen, "&"))
 # Sys.sleep(2)
#} # if not strand extract
#






######################################################
################## PREPARE GRAPHS ####################
######################################################
setwd(collectDir)
featureCounts.counts <- read.csv(featureCounts.counts.fName, header = TRUE, row.names=1)
featureCounts.FPKM <- read.csv(featureCounts.RPKM.fName, header = TRUE, row.names=1)
# Create graphs
featureCounts.counts.graphTitle <- paste(text.add,",log10(counts+1",sep="")
featureCounts.FPKM.graphTitle <- paste(text.add,",log10(FPKM+0.01",sep="")
#
graph.title <- paste(text.add,".FeatureCounts.BoxPlot.png",sep="")

png(graph.title,width=1400,height=700)
par(mfrow=c(1,2))
boxplot(log10(featureCounts.counts + 1), main=featureCounts.counts.graphTitle)
boxplot(log10(featureCounts.FPKM + 0.01), main=featureCounts.FPKM.graphTitle)
dev.off()