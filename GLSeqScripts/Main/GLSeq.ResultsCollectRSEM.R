args <- commandArgs(trailingOnly = TRUE)
load.file <- as.character(args[1])
load(load.file)


setwd(destDirRSEMCount)
bai.pull <- function(dest.dir, text.add) {
  allfiles <- dir(dest.dir)
  bai.files <- allfiles[grep("bam.bai", allfiles)]
  # just in case, let's restrict the list to the files with the current "text.add" mark:
  bai.files.thisRun <- bai.files[grep(text.add, bai.files)]
  # those are mix of transcript- and genome-level bai files;
  # restricting to genome-level files:
  bai.files.thisRun.genome <- bai.files.thisRun[grep("genome.sorted.bam.bai", bai.files.thisRun)]
  # if genome bam file is not requested, then transcriptome bai (the only indexed bam's captured by bai.pull() )is being used instead:
  if (!(genobam)) bai.files.thisRun.genome <- bai.files
  bai.files.thisRun.genome
}
#
############################
# generating matrices of computation results
# for all the libraries in the run
############################
#
# name signature of the results files:
results.fileSig <- paste("genes.results", sep=".")
# file names with library-centric results:
result.fnames <- dir(pattern=results.fileSig)
print("RESULTS")
print(result.fnames)
print(length(result.fnames))
# library name + text add:
result.names <- substr(result.fnames, 1, nchar(result.fnames)-14)
# just library names, in the respective order:
lib.names <- substr(result.names, 1,libNchar)
# collecting expected count and other processed data:
counts <- NULL
counts_pme <- NULL
FPKM <- NULL
FPKM_pme <- NULL
FPKM_lower <- NULL
FPKM_upper <- NULL
TPM <- NULL
TPM_pme <- NULL
TPM_lower <- NULL
TPM_upper <- NULL
#
for (i in 1:length(result.fnames)) {
  i.data <- read.table(result.fnames[i], header=TRUE, sep="\t", row.names=1, as.is=TRUE)
  counts <- cbind(counts, i.data[,"expected_count"])
  # counts_pme <- cbind(counts_pme, i.data[,"pme_expected_count"]) # field name changed in RSEM v.1.2.15 (Jun 16, 2014)
  counts_pme <- cbind(counts_pme, i.data[,"posterior_mean_count"])
  FPKM <- cbind(FPKM, i.data[,"FPKM"])
  FPKM_pme <- cbind(FPKM_pme, i.data[,"pme_FPKM"])
  FPKM_lower <- cbind(FPKM_lower, i.data[,"FPKM_ci_lower_bound"])
  FPKM_upper <- cbind(FPKM_upper, i.data[,"FPKM_ci_upper_bound"])
  TPM <- cbind(TPM, i.data[,"TPM"])
  TPM_pme <-  cbind(TPM_pme, i.data[,"pme_TPM"])
  TPM_lower <-  cbind(TPM_lower, i.data[,"TPM_ci_lower_bound"])
  TPM_upper <- cbind(TPM_upper, i.data[,"TPM_ci_upper_bound"])
} # for i
#
colnames(counts) <- result.names
colnames(counts_pme) <- result.names
colnames(FPKM) <- result.names
colnames(FPKM_pme) <- result.names
colnames(FPKM_lower) <- result.names
colnames(FPKM_upper) <- result.names
colnames(TPM) <- result.names
colnames(TPM_pme) <- result.names
colnames(TPM_lower) <- result.names
colnames(TPM_upper) <- result.names
#
rownames(counts) <- rownames(i.data)
rownames(counts_pme) <- rownames(i.data)
rownames(FPKM) <- rownames(i.data)
rownames(FPKM_pme) <- rownames(i.data)
rownames(FPKM_lower) <- rownames(i.data)
rownames(FPKM_upper) <- rownames(i.data)
rownames(TPM) <- rownames(i.data)
rownames(TPM_pme) <- rownames(i.data)
rownames(TPM_lower) <- rownames(i.data)
rownames(TPM_upper) <- rownames(i.data)
#
counts <- round(counts, 0)
counts_pme <- round(counts_pme, 0)
#
counts.fName <- paste(collectDir, text.add, ".counts.csv", sep="")
counts_pme.fName <-  paste(collectDir, text.add, ".counts_pme.csv", sep="")
FPKM.fName <-  paste(collectDir, text.add, ".FPKM.csv", sep="")
FPKM_pme.fName <-  paste(collectDir, text.add, ".FPKM_pme.csv", sep="")
FPKM_lower.fName <-  paste(collectDir, text.add, ".FPKM_lower.csv", sep="")
FPKM_upper.fName <-  paste(collectDir, text.add, ".FPKM_upper.csv", sep="")
TPM.fName <-  paste(collectDir, text.add, ".TPM.csv", sep="")
TPM_pme.fName <-  paste(collectDir, text.add, ".TPM_pme.csv", sep="")
TPM_lower.fName <-  paste(collectDir, text.add, ".TPM_lower.csv", sep="")
TPM_upper.fName <-  paste(collectDir, text.add, ".TPM_upper.csv", sep="")
#
write.csv(counts, file=counts.fName)
write.csv(counts_pme, file=counts_pme.fName)
write.csv(FPKM, file=FPKM.fName)
write.csv(FPKM_pme, file=FPKM_pme.fName)
write.csv(FPKM_lower, file=FPKM_lower.fName)
write.csv(FPKM_upper, file=FPKM_upper.fName)
write.csv(TPM, file=TPM.fName)
write.csv(TPM_pme, file=TPM_pme.fName)
write.csv(TPM_lower, file=TPM_lower.fName)
write.csv(TPM_upper, file=TPM_upper.fName)
#
############################
# collecting genome-level BAM and BAI files,
# generating wiggle files
# and moving all of them (.bam, .bai, .wig)
# to a chosen folder (destDirBam)
############################
# if genome file is not requested (genobam == FALSE), then all the files below will indeed represent reanscriptome visualization files
# because bai.files.thisRun.genome will actually contain transcriptome bam's (see bai.pull function for details)
bai.files.thisRun.genome <- bai.pull(dest.dir, text.add)
# respective names of the bam and wig files:
bam.files.thisRun.genome <- substr(bai.files.thisRun.genome, 1, nchar(bai.files.thisRun.genome)-4)
wig.files.thisRun.genome <- paste(substr(bam.files.thisRun.genome, 1, nchar(bam.files.thisRun.genome)-4),"wig", sep=".")
#
# The case for generating wigle files without explicit separating forward- and reverse-strand reads
# (implies using of a standard gtf file without addition of reverse-strand features):
if(!(strandExtract)) {
  # Generating .wig files:
  for (bam in 1:length(bam.files.thisRun.genome)) {
    current.bam <- bam.files.thisRun.genome[bam]
    current.bai <- bai.files.thisRun.genome[bam]
    # name of the wig file to generate:
    current.wig <- wig.files.thisRun.genome[bam]
    # base name (before ".genome.sorted.bam"):
    current.base <- substr(current.bam, 1, nchar(current.bam)-18)
    wig.command <- paste("rsem-bam2wig ", current.bam, current.wig, current.base, sep=" ")
    system(wig.command)
  }
  cat(timestamp(), " Wiggle files are generated!", "\n")
} # if not strandExtract
#
######################################################
################## PREPARE GRAPHS ####################
######################################################
setwd(collectDir)
RSEM.counts <- read.csv(counts_pme.fName, header = TRUE, row.names=1)
RSEM.FPKM <- read.csv(FPKM_pme.fName, header = TRUE, row.names=1)
#
RSEM.counts.graphTitle <- "RSEM log10(Counts+1)"
RSEM.FPKM.graphTitle <- "RSEM log10(FPKM+0.01)"
#
graph.title <- paste(text.add,".RSEM.BoxPlot.png",sep="")
png(graph.title,width=1400,height=700)
par(mfrow=c(1,2))
boxplot(log10(RSEM.counts + 1), main=RSEM.counts.graphTitle)
boxplot(log10(RSEM.FPKM + 0.01), main=RSEM.FPKM.graphTitle)
dev.off()