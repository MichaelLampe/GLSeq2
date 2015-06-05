#########################################################
# Great Lakes Seq package for low-level processing of RNA-Seq data
#########################################################
#
# This is a user generated attribute file created on 2015.06.05
#
##########################################################
#
# Usage: called from the GLSeq.top.R script;
# multiple versions of local attribute files may exist; the name of the particular version
# may be supplied as an option to the GLSeq.top.R script
#
##########################################################
#
################################
# DATA / LIBRARY OPTIONS
###############################
#
# directory containing raw files
# (may be non-writable!)
raw.dir <-""
#
# Files in the raw dir are normally compressed but may be not:
unzipped <- FALSE
#
# directory contining ready-to-go (split+QC-processed) fq files (Oct 17, 2013)
readyData.dir <- "/home/GLBRCORG/mrlampe/GLBRC_UI/Testing/.Test_Cases/Test_CushawGPU_527_HTSeq_FeatureCounts.05"
#
# raw file names: 
raw.fNames <- ""
#
# strain
strain <- ""
#
# single / paired end
paired.end <- TRUE
#
# sequencing platform (used by CUSHAW, supported values: capillary, ls454, illumina, solid, helicos, iontorrent, pacbio)
seqPlatform <- "illumina"
#
# quality scores format
qScores <- "phred33"
#
# Strandness of the library (NULL, F, R) 
libstrand <- "R"
#
# Number of unique characters in the beginning of the each file (library ID length):
libNchar <- 4
#
# Subset of the libraries to process (optional; normally the list wil be generated from the actual directory content)
libList <- "NULL"
#
# Takes a directory of files with the end title "countable.sam" and collects them for counting
countable.sams.dir <- ""
#
###############################
# REFERENCE OPTIONS
###############################
#
# reference genome - the index directory for the respective method (RSEM or BWA)
# (must match the name of the  subfolder under base.dir):
rGenome <- "Y22-3.A144"
#
# name of the reference fasta file (may differ from the base name of the reference -
# still, it should be located in the folder for the selected feference):
refFASTAname <- "Y22-3.assembly.2015_01_13.fasta"
#
# name of the reference genomic features file (may differ from the base name of the reference -
# this will be eventually taken from the database); it should be located in the folder for the selected feference)
# gtf files may be used instead of gff where applicable; the same object is used for both cases:
refGFFname <- "Y22-3.A144c.gtf"
# refGFFname <- "Novosphingobium_3replicons.Clean.gff"
#
# number of the column in GTF file with the gene / other IDs charachter string (9, unless the file is non-standard for some reason):
gtfFeatureColumn <- 9
#
# GFF attribute to be used as feature ID (HTSeq):
idAttr <- "gene_id"
# idAttr <- "locus_tag" # useful if "gene_id" has duplicated entries
#
###############################
# RUN OPTIONS
###############################
#
# the directory containing the GLSeq scripts:
base.dir <- "/home/GLBRCORG/mrlampe/GLBRC_UI/Working_RScript_Files"
#
# Base of the destination directory (added May 9, 2013)
# This should be located on a FAST volume (SCSI is recommended)
# a particular subfolder names after the run ID will be created by GLSeq below this folder
dest.dir.base <- "/home/GLBRCORG/mrlampe/GLBRC_UI/Testing"
#
# number of cores to use
nCores <- 7
#
# number of parallel computation streams for expression computation
nStreams <-3
#
# number of parallel computation streams for data preparation
# (may differ from the number of streams for expression computation because of particular software demands) 
nStreamsDataPrep <- 1
#
# the actual unique run ID - 
# text.add <- paste(expID, runAttempt, sep=".")
# now is being generated inside GLSeq.top.R
#
# *** Alignment Algorithm ***
aAlgor <- "Cushaw" 
#
#
# *** Counting Algorithm(s) ***
HTSeq <- "HTSeq"
FeatureCounts <- "FeatureCounts"
RSEM <- ""
cAlgor <- c(HTSeq,RSEM,FeatureCounts)
#
#  GPU acceleration option for CUSHAW
GPU.accel <- TRUE
#
###############################
# PRE-PROCESSING OPTIONS
###############################
#
# trim the reads and generate QC reports for before- and after-trimming FASTQ files? 
readTrim <- TRUE
#
# minimum length of a trimmed read
trimMin <- 0
#
# trimmomatic parameter values for HEADCROP  
trimhead <- 0
#
# name of the FASTA file with artificail sequences (adapters, primers etc) - must be located in the base.dir
artificial.fq <- ""
#
###############################
# COMMON PROCESSING OPTIONS
###############################
#
# Extract explicit forward and reverse coverage from the original BAM file? 
strandExtract <- FALSE
#
################################
# RSEM OPTIONS
################################
#
# compute confidence intervals? 
compConf <- FALSE
#
# Maximal length of fragment (for paired-end libraries)
fragMaxLength <- 1000
#
# Maximum size (MB) of the auxiliary buffer used for computing credibility intervals (CI) - for RSEM (+extra 2Gb per stream)
ciMem <- 4096
#
# Output genome bam
genobam <- FALSE
#
################################
# ENVIRONMENT
################################
#
# path to Trimmomatic: 
trimPath <- '/opt/bifxapps/Trimmomatic-0.30/trimmomatic-0.30.jar'
#
# path to PicardTools jar directory:
# picardToolsPath <- '/soft/picard-tools-1.98/picard-tools-1.98/'
picardToolsPath <- '/opt/bifxapps/picard-tools/'
#
# path to fastqc:
fastqcPath <- '/opt/bifxapps/bin/fastqc'
#
# path to BWA
bwaPath <- '/opt/bifxapps/bin/bwa'
#
# path to the shell script that converts bam to wig
bam2wigPath <- "/home/GLBRCORG/omoskvin/run/bam2wig.sh"
#
# path to CUSHAW
CUSHAW.index.path <- "/opt/bifxapps/cushaw2-v2.1.11/cushaw2_index/cushaw2_index"
CUSHAW.path <- "/opt/bifxapps/cushaw2-v2.1.11/cushaw2"
#
# path to CUSHAW-GPU
CUSHAW.GPU.path <- "/opt/bifxapps/cushaw2-gpu-2.1.8-r16/cushaw2-gpu"
#
# End of Attribute File
