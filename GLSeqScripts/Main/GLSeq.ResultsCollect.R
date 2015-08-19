############
# GLSeq2 Collection
############
#
setwd(dest.dir)
# Otherwise the preloaded ones are fine.
#
#
#







































#############################################################################################
###############################             HT SEQ            ###############################
#############################################################################################
if ("HTSeq" %in% cAlgor){
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
}
#
#





























#############################################################################################
###############################         FEATURE COUNTS        ###############################
#############################################################################################
if ("FeatureCounts" %in% cAlgor){
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
  cfiles <- dir(pattern="counts")
  countDirName <- paste(text.add, "counts", sep=".")
  cfiles <- cfiles[cfiles != countDirName]
  count.mtrx <- NULL
  #
  #
  # Need to fix the formatting of the file because FeatureCounts has a bunch of useless/artistic outputs that clutter the file.
  #
  correctedNames <- c()
  for (i in 1:length(cfiles)){
    setwd(destDirFeatureCountsCount)
    correctedName = paste(cfiles[i],".corrected",sep="")
    correctedNames <- c(correctedNames,correctedName)
    openFile <- file(cfiles[i], open="r")
    document <- readLines(openFile)
    for (i in 1:length(document))
      if (grepl("$counts",document[i],fixed=TRUE)){
        i <- i + 2
        setwd(collectDir)
        while (!(grepl("$annotation",document[i],fixed=TRUE))){
          write(document[i],file=correctedName,append=TRUE)
          i <- i + 1
        }
        currentLine <- readLines(openFile)
      }
    close(openFile)
  }
  cfiles <- correctedNames
  #
  for (i in 1:length(cfiles)) {
    data.i <- read.table(cfiles[i], row.names=1)
    colnames(data.i) <- substr(cfiles[i],1,libNchar)
    if (is.null(count.mtrx)) count.mtrx <- data.i
    if (!(is.null(count.mtrx)) & i != 1) count.mtrx <- cbind(count.mtrx, data.i)
  }
  #
  featureCounts.counts.fName <-  paste(collectDir, text.add, ".FeatureCounts.counts.csv", sep="")
  write.csv(count.mtrx, file=featureCounts.counts.fName)
  #
  #
  # Load GTF file
  #
  reference.gtf <- paste(dest.dir,refGFFname,sep="")
  gtf <-  read.table(reference.gtf, sep="\t", header=FALSE, as.is=TRUE)
  mess <- gtf[,gtfFeatureColumn] # 9
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
  #
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
  } # if not strand extract
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
} # Feature Counts
#
#
#
















































#############################################################################################
###############################            RSEM               ###############################
#############################################################################################
if ("RSEM" %in% cAlgor){
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
} # RSEM