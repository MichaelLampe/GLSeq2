package org.glbrc.glseq2;
///////////////////////////////////////////////////////////////////////////////
//                   
// Main Class File:  UI.java
// File:             Attributes.java  
//
//
// Author:           Michael Lampe | MrLampe@Wisc.edu
///////////////////////////////////////////////////////////////////////////////
import java.io.*;
import java.lang.reflect.*;
import java.util.*;
import java.text.*;

public class Attributes {

	// Data and Library Attributes
	private String directory = "";
	private String unzipped = "FALSE";
	private String directoryFq = "";
	private String rawFileNames = "";
	private String strain = "";
	private String pairedEnd = "TRUE";
	private String seqPlatform = "";
	private String qScores = "";
	private String libstrands = "";
	private String libNchar = "0";
	private String libList = "\"NULL\"";
	private String countableSamDir = "";

	// Reference Attributes
	private String rGenome = "";
	private String refFASTAname = "";
	private String refGFFname = "";
	private String gtfFeatureColumn = "0";
	private String idAttr = "";

	// Run Attributes
	private String destinationDirectory = "";
	private String FeatureCounts = "";
	private String RSEM = "";
	private String HTSeq = "";
	private String aAlgor = "";
	private String nCores = "0";
	private String nStreams = "0";
	private String nStreamsDataPrep = "0";
	private String GPUaccel = "FALSE";

	// Pre-Processing Attributes
	private String readTrim = "0";
	private String trimHead = "0";
	private String minTrim = "0";
	private String artificialFASTA = "";

	// Common Processing Attributes
	private String strandExtract = "FALSE";

	// RSEM Attributes
	private String compConf = "FALSE";
	private String fragMaxLength = "0";
	private String ciMem = "0";
	private String genobam = "FALSE";

	// Environment Attributes
	private String scriptDirectory = "";
	private String trimPath = "";
	private String picardToolsPath = "";
	private String fastqcPath = "";
	private String bwaPath = "";
	private String bam2wigPath = "";
	private String cushawPath = "";
	private String cushawGpuPath = "";
	private String cushawIndexPath = "";

	// RunId
	public String runId = "";

	/**
	 * Creates a new or replaces the current attribute file and adds all the
	 * field data to that attribute file.
	 * 
	 * @throws IllegalAccessException
	 * @throws IllegalArgumentException
	 * @throws SecurityException
	 * @throws NoSuchFieldException
	 *
	 */
	public void setAttributes() throws SecurityException {
		Field field = null;
		Scanner scnr = null;
		File config = null;
		try {
			config = new File("AttributesConfig.txt");
			scnr = new Scanner(config);
		} catch (Exception e) {
		}
		// Looks through the whole document for variables that match field names
		while (scnr.hasNextLine()) {
			String temp = scnr.nextLine();
			if (temp.contains("= ")) {
				String[] t = temp.split(" = ");
				try {
				field = this.getClass().getDeclaredField(t[0]);
				} catch (NoSuchFieldException e){
				}
				try {
					field.set(this, (t[1]));
				} catch (IllegalAccessException e) {
					e.printStackTrace();
				} catch (IllegalArgumentException e){
					e.printStackTrace();
				} catch(ArrayIndexOutOfBoundsException e){
				}
			}
		}
		scnr.close();
	}

	/**
	 * Creates a new or replaces the current attribute file and adds all the
	 * field data to that attribute file.
	 * 
	 * This method uses an input attributeFile.R to take old data and load it
	 * in.
	 *
	 */
	public void setAttributes(File attributeFile) {
		Field field;
		Scanner scnr = null;
		try {
			scnr = new Scanner(attributeFile);
			// Looks through the whole document for variables that match field
			// names
			while (scnr.hasNextLine()) {
				String temp = scnr.nextLine();
				if (temp.contains("=")) {
					String[] t = temp.split(" = ");
					try {
						field = this.getClass().getDeclaredField(t[0]);
						field.set(this, t[1]);
					} catch (Exception e) {
					}
				}
			}
		} catch (Exception e) {
			System.out.println("Attribute file not found");
		} finally {
			// Makes sure to close the scanner
			scnr.close();
		}
	}

	// Overloaded for user to provide file location
	public void setAttributes(String s) {
	}

	// Setters for all the variables.
	public void setDirectory(String directory) {
		this.directory = directory;
	}

	public void setUnzipped(String unzipped) {
		this.unzipped = unzipped;
	}

	public void setDirectoryFq(String directoryFq) {
		this.directoryFq = directoryFq;
	}

	public void setRawFileNames(String rawFileNames) {
		this.rawFileNames = rawFileNames;
	}

	public void setStrain(String strain) {
		this.strain = strain;
	}

	public void setPairedEnd(String pairedEnd) {
		this.pairedEnd = pairedEnd;
	}

	public void setSeqPlatform(String seqPlatform) {
		this.seqPlatform = seqPlatform;
	}

	public void setQScores(String qScores) {
		this.qScores = qScores;
	}

	public void setLibstrand(String libstrands) {
		this.libstrands = libstrands;
	}

	public void setLibNchar(String libNchar) {
		this.libNchar = libNchar;
	}

	public void setLibList(String libList) {
		this.libList = libList;
	}
	public void setCountableSamDir(String countableSamDir){
		this.countableSamDir = countableSamDir;
	}

	// Reference Attribute Setters
	public void setRGenome(String rGenome) {
		this.rGenome = rGenome;
	}

	public void setFASTAname(String refFASTAname) {
		this.refFASTAname = refFASTAname;
	}

	public void setGFFname(String refGFFname) {
		this.refGFFname = refGFFname;
	}

	public void setFeatureColumn(String gtfFeatureColumn) {
		this.gtfFeatureColumn = gtfFeatureColumn;
	}

	public void setIdAttr(String idAttr) {
		this.idAttr = idAttr;
	}

	// Run Attributes
	public void setScriptDirectory(String scriptDirectory) {
		this.scriptDirectory = scriptDirectory;
	}

	public void setDestinationDirectory(String destinationDirectory) {
		this.destinationDirectory = destinationDirectory;
	}

	public void setCores(String nCores) {
		this.nCores = nCores;
	}

	public void setStreams(String nStreams) {
		this.nStreams = nStreams;
	}

	public void setStreamsDataPrep(String nStreamsDataPrep) {
		this.nStreamsDataPrep = nStreamsDataPrep;
	}

	public void setGpuAccel(String GPUaccel) {
		this.GPUaccel = GPUaccel;
	}

	// Pre-Processing Attributes
	public void setReadTrim(String readTrim) {
		this.readTrim = readTrim;
	}

	public void setTrimhead(String trimHead) {
		this.trimHead = trimHead;
	}

	public void setMinTrim(String minTrim) {
		this.minTrim = minTrim;
	}

	public void setArtificialFASTA(String artificialFASTA) {
		this.artificialFASTA = artificialFASTA;
	}

	// Common Processing Attributes
	public void setStrandExtract(String strandExtract) {
		this.strandExtract = strandExtract;
	}

	// RSEM Attributes
	public void setCompConf(String compConf) {
		this.compConf = compConf;
	}

	public void setFragMaxLength(String fragMaxLength) {
		this.fragMaxLength = fragMaxLength;
	}

	public void setCiMem(String ciMem) {
		this.ciMem = ciMem;
	}

	public void setGenoBam(String genobam) {
		this.genobam = genobam;
	}

	// Environment Attributes
	public void setTrimPath(String trimPath) {
		this.trimPath = trimPath;
	}

	public void setPicardToolsPath(String picardToolsPath) {
		this.picardToolsPath = picardToolsPath;
	}

	public void setFastqcPath(String fastqcPath) {
		this.fastqcPath = fastqcPath;
	}

	public void setBwaPath(String bwaPath) {
		this.bwaPath = bwaPath;
	}

	public void setBam2WigPath(String bam2wigPath) {
		this.bam2wigPath = bam2wigPath;
	}

	public void setCushawPath(String cushawPath) {
		this.cushawPath = cushawPath;
	}

	public void setCushawGpuPath(String cushawGpuPath) {
		this.cushawGpuPath = cushawGpuPath;
	}
	//
	//
	public void setCushawIndexPath(String cushawIndexPath){
		this.cushawIndexPath = cushawIndexPath;
	}
	public void setHTSeq(String HTSeq){
		this.HTSeq= HTSeq;
	}
	public void setFeatureCounts(String FeatureCounts){
		this.FeatureCounts = FeatureCounts;
	}
	public void setRSEM(String RSEM){
		this.RSEM = RSEM;
	}
	public void setaAlgor(String aAlgor){
		this.aAlgor = aAlgor;
	}
	
	//
	public String getCushawIndexPath(){
		return cushawIndexPath;
	}
	public String getHTSeq(){
		return HTSeq;
	}
	public String getFeatureCounts(){
		return FeatureCounts;
	}
	public String getRSEM(){
		return RSEM;
	}
	public String getaAlgor(){
		return aAlgor;
	}

	// Data and Library Getters
	public String getDirectory() {
		return directory;
	}

	public String getUnzipped() {
		return unzipped;
	}

	public String getDirectoryFq() {
		return directoryFq;
	}

	public String getRawFileNames() {
		return rawFileNames;
	}

	public String getStrain() {
		return strain;
	}

	public String getPairedEnd() {
		return pairedEnd;
	}

	public String getSeqPlatform() {
		return seqPlatform;
	}

	public String getQScores() {
		return qScores;
	}

	public String getLibstrand() {
		return libstrands;
	}

	public String getLibNchar() {
		return libNchar;
	}

	public String getLibList() {
		return libList;
	}
	public String getCountableSamDir(){
		return countableSamDir;
	}

	// Reference Attribute Getters
	public String getRGenome() {
		return rGenome;
	}

	public String getFASTAname() {
		return refFASTAname;
	}

	public String getGFFname() {
		return refGFFname;
	}

	public String getFeatureColumn() {
		return gtfFeatureColumn;
	}

	public String getIdAttr() {
		return idAttr;
	}

	// Run Attributes
	public String getScriptDirectory() {
		return scriptDirectory;
	}

	public String getDestinationDirectory() {
		return destinationDirectory;
	}

	public String getCores() {
		return nCores;
	}

	public String getStreams() {
		return nStreams;
	}

	public String getStreamsDataPrep() {
		return nStreamsDataPrep;
	}

	public String getGpuAccel() {
		return GPUaccel;
	}
	// Pre-Processing Attributes
	public String getReadTrim() {
		return readTrim;
	}

	public String getTrimhead() {
		return trimHead;
	}

	public String getMinTrim() {
		return minTrim;
	}

	public String getArtificialFASTA() {
		return artificialFASTA;
	}

	// Common Processing Attributes
	public String getStrandExtract() {
		return strandExtract;
	}

	// RSEM Attributes
	public String getCompConf() {
		return compConf;
	}

	public String getFragMaxLength() {
		return fragMaxLength;
	}

	public String getCiMem() {
		return ciMem;
	}

	public String getGenoBam() {
		return genobam;
	}

	// Environment Attributes
	public String getTrimPath() {
		return trimPath;
	}

	public String getPicardToolsPath() {
		return picardToolsPath;
	}

	public String getFastqcPath() {
		return fastqcPath;
	}

	public String getBwaPath() {
		return bwaPath;
	}

	public String getBam2WigPath() {
		return bam2wigPath;
	}

	public String getCushawPath() {
		return cushawPath;
	}

	public String getCushawGpuPath() {
		return cushawGpuPath;
	}

	public String writeAttributesFile() throws IOException {
		// Writes a new attribute file based on the current format that we have
		// outlined.
		// Saves it as an R file, with the current date.
		if (aAlgor.equals("Cushaw-GPU")){
			aAlgor = "Cushaw";
		}
		Date current = new Date();
		SimpleDateFormat ft = new SimpleDateFormat("yyyy.MM.dd");
		File attributeFile = new File("Attribute_" + ft.format(current) + ".R");
		attributeFile.createNewFile();
		String attFileLocation = attributeFile.getAbsolutePath();
		FileWriter writer = new FileWriter(attributeFile);
		writer.write("#########################################################\n");
		writer.write("# Great Lakes Seq package for low-level processing of RNA-Seq data\n");
		writer.write("#########################################################\n");
		writer.write("#\n");
		writer.write("# This is a user generated attribute file created on "
				+ ft.format(current) + "\n");
		writer.write("#\n");
		writer.write("##########################################################\n");
		writer.write("#\n");
		writer.write("# Usage: called from the GLSeq.top.R script;\n");
		writer.write("# multiple versions of local attribute files may exist; the name of the particular version\n");
		writer.write("# may be supplied as an option to the GLSeq.top.R script\n");
		writer.write("#\n");
		writer.write("##########################################################\n");
		writer.write("#\n");
		writer.write("################################\n");
		writer.write("# DATA / LIBRARY OPTIONS\n");
		writer.write("###############################\n");
		writer.write("#\n");
		writer.write("# directory containing raw files\n");
		writer.write("# (may be non-writable!)\n");
		writer.write("raw.dir <-\"" + directory + "\"\n");
		writer.write("#\n");
		writer.write("# Files in the raw dir are normally compressed but may be not:\n");
		writer.write("unzipped <- " + unzipped.toUpperCase() + "\n");
		writer.write("#\n");
		writer.write("# directory contining ready-to-go (split+QC-processed) fq files (Oct 17, 2013)\n");
		writer.write("readyData.dir <- \"" + directoryFq + "\"\n");
		writer.write("#\n");
		writer.write("# raw file names: \n");
		writer.write("raw.fNames <- \"" + rawFileNames + "\"\n");
		writer.write("#\n");
		writer.write("# strain\n");
		writer.write("strain <- \"" + strain + "\"\n");
		writer.write("#\n");
		writer.write("# single / paired end\n");
		writer.write("paired.end <- " + pairedEnd.toUpperCase() + "\n");
		writer.write("#\n");
		writer.write("# sequencing platform (used by CUSHAW, supported values: capillary, ls454, illumina, solid, helicos, iontorrent, pacbio)\n");
		writer.write("seqPlatform <- \"" + seqPlatform + "\"\n");
		writer.write("#\n");
		writer.write("# quality scores format\n");
		writer.write("qScores <- \"" + qScores + "\"\n");
		writer.write("#\n");
		writer.write("# Strandness of the library (NULL, F, R) \n");
		writer.write("libstrand <- \"" + libstrands + "\"\n");
		writer.write("#\n");
		writer.write("# Number of unique characters in the beginning of the each file (library ID length):\n");
		writer.write("libNchar <- " + libNchar + "\n");
		writer.write("#\n");
		writer.write("# Subset of the libraries to process (optional; normally the list wil be generated from the actual directory content)\n");
		writer.write("libList <- " + libList + "\n");
		writer.write("#\n");
		writer.write("# Takes a directory of files with the end title \"countable.sam\" and collects them for counting\n");
		writer.write("countable.sams.dir <- \"" + countableSamDir + "\"\n");
		writer.write("#\n");
		writer.write("###############################\n");
		writer.write("# REFERENCE OPTIONS\n");
		writer.write("###############################\n");
		writer.write("#\n");
		writer.write("# reference genome - the index directory for the respective method (RSEM or BWA)\n");
		writer.write("# (must match the name of the  subfolder under base.dir):\n");
		writer.write("rGenome <- \"" + rGenome + "\"\n");
		writer.write("#\n");
		writer.write("# name of the reference fasta file (may differ from the base name of the reference -\n");
		writer.write("# still, it should be located in the folder for the selected feference):\n");
		writer.write("refFASTAname <- \"" + refFASTAname + "\"\n");
		writer.write("#\n");
		writer.write("# name of the reference genomic features file (may differ from the base name of the reference -\n");
		writer.write("# this will be eventually taken from the database); it should be located in the folder for the selected feference)\n");
		writer.write("# gtf files may be used instead of gff where applicable; the same object is used for both cases:\n");
		writer.write("refGFFname <- \"" + refGFFname + "\"\n");
		writer.write("# refGFFname <- \"Novosphingobium_3replicons.Clean.gff\"\n");
		writer.write("#\n");
		writer.write("# number of the column in GTF file with the gene / other IDs charachter string (9, unless the file is non-standard for some reason):\n");
		writer.write("gtfFeatureColumn <- " + gtfFeatureColumn + "\n");
		writer.write("#\n");
		writer.write("# GFF attribute to be used as feature ID (HTSeq):\n");
		writer.write("idAttr <- \"" + idAttr + "\"\n");
		writer.write("# idAttr <- \"locus_tag\" # useful if \"gene_id\" has duplicated entries\n");
		writer.write("#\n");
		writer.write("###############################\n");
		writer.write("# RUN OPTIONS\n");
		writer.write("###############################\n");
		writer.write("#\n");
		writer.write("# the directory containing the GLSeq scripts:\n");
		writer.write("base.dir <- \"" + scriptDirectory + "\"\n");
		writer.write("#\n");
		writer.write("# Base of the destination directory (added May 9, 2013)\n");
		writer.write("# This should be located on a FAST volume (SCSI is recommended)\n");
		writer.write("# a particular subfolder names after the run ID will be created by GLSeq below this folder\n");
		writer.write("dest.dir.base <- \"" + destinationDirectory + "\"\n");
		writer.write("#\n");
		writer.write("# number of cores to use\n");
		writer.write("nCores <- " + nCores + "\n");
		writer.write("#\n");
		writer.write("# number of parallel computation streams for expression computation\n");
		writer.write("nStreams <-" + nStreams + "\n");
		writer.write("#\n");
		writer.write("# number of parallel computation streams for data preparation\n");
		writer.write("# (may differ from the number of streams for expression computation because of particular software demands) \n");
		writer.write("nStreamsDataPrep <- " + nStreamsDataPrep + "\n");
		writer.write("#\n");
		writer.write("# the actual unique run ID - \n");
		writer.write("# text.add <- paste(expID, runAttempt, sep=\".\")\n");
		writer.write("# now is being generated inside GLSeq.top.R\n");
		writer.write("#\n");
		writer.write("# *** Alignment Algorithm ***\n");
		writer.write("aAlgor <- \"" + aAlgor + "\" \n");
		writer.write("#\n");
		writer.write("#\n");
		// This is just assigning them as variables and just making them as a list in the R file
		writer.write("# *** Counting Algorithm(s) ***\n");
		writer.write("HTSeq <- \"" +HTSeq + "\"\n");
		writer.write("FeatureCounts <- \"" + FeatureCounts + "\"\n");
		writer.write("RSEM <- \"" + RSEM + "\"\n");
		writer.write("cAlgor <- c(HTSeq,RSEM,FeatureCounts)\n");
		writer.write("#\n");
		writer.write("#  GPU acceleration option for CUSHAW\n");
		writer.write("GPU.accel <- " + GPUaccel.toUpperCase() + "\n");
		writer.write("#\n");
		writer.write("###############################\n");
		writer.write("# PRE-PROCESSING OPTIONS\n");
		writer.write("###############################\n");
		writer.write("#\n");
		writer.write("# trim the reads and generate QC reports for before- and after-trimming FASTQ files? \n");
		writer.write("readTrim <- " + readTrim.toUpperCase() + "\n");
		writer.write("#\n");
		writer.write("# minimum length of a trimmed read\n");
		writer.write("trimMin <- " + minTrim + "\n");
		writer.write("#\n");
		writer.write("# trimmomatic parameter values for HEADCROP  \n");
		writer.write("trimhead <- " + trimHead + "\n");
		writer.write("#\n");
		writer.write("# name of the FASTA file with artificail sequences (adapters, primers etc) - must be located in the base.dir\n");
		writer.write("artificial.fq <- \"" + artificialFASTA + "\"\n");
		writer.write("#\n");
		writer.write("###############################\n");
		writer.write("# COMMON PROCESSING OPTIONS\n");
		writer.write("###############################\n");
		writer.write("#\n");
		writer.write("# Extract explicit forward and reverse coverage from the original BAM file? \n");
		writer.write("strandExtract <- " + strandExtract.toUpperCase() + "\n");
		writer.write("#\n");
		writer.write("################################\n");
		writer.write("# RSEM OPTIONS\n");
		writer.write("################################\n");
		writer.write("#\n");
		writer.write("# compute confidence intervals? \n");
		writer.write("compConf <- " + compConf.toUpperCase() + "\n");
		writer.write("#\n");
		writer.write("# Maximal length of fragment (for paired-end libraries)\n");
		writer.write("fragMaxLength <- " + fragMaxLength + "\n");
		writer.write("#\n");
		writer.write("# Maximum size (MB) of the auxiliary buffer used for computing credibility intervals (CI) - for RSEM (+extra 2Gb per stream)\n");
		writer.write("ciMem <- " + ciMem + "\n");
		writer.write("#\n");
		writer.write("# Output genome bam\n");
		writer.write("genobam <- " + genobam.toUpperCase() + "\n");
		writer.write("#\n");
		writer.write("################################\n");
		writer.write("# ENVIRONMENT\n");
		writer.write("################################\n");
		writer.write("#\n");
		writer.write("# path to Trimmomatic: \n");
		writer.write("trimPath <- '" + trimPath + "'\n");
		writer.write("#\n");
		writer.write("# path to PicardTools jar directory:\n");
		writer.write("# picardToolsPath <- '/soft/picard-tools-1.98/picard-tools-1.98/'\n");
		writer.write("picardToolsPath <- '" + picardToolsPath + "'\n");
		writer.write("#\n");
		writer.write("# path to fastqc:\n");
		writer.write("fastqcPath <- '" + fastqcPath + "'\n");
		writer.write("#\n");
		writer.write("# path to BWA\n");
		writer.write("bwaPath <- '" + bwaPath + "'\n");
		writer.write("#\n");
		writer.write("# path to the shell script that converts bam to wig\n");
		writer.write("bam2wigPath <- \"" + bam2wigPath + "\"\n");
		writer.write("#\n");
		writer.write("# path to CUSHAW\n");
		writer.write("CUSHAW.index.path <- \"" + cushawIndexPath + "\"\n");
		writer.write("CUSHAW.path <- \"" + cushawPath + "\"\n");
		writer.write("#\n");
		writer.write("# path to CUSHAW-GPU\n");
		writer.write("CUSHAW.GPU.path <- \"" + cushawGpuPath + "\"\n");
		writer.write("#\n");
		writer.write("# End of Attribute File\n");
		writer.close();
		return attFileLocation;
	}

	public void saveConfigFile(String fileName) throws IOException {
		// An array of all the fields
		Field[] field = this.getClass().getDeclaredFields();
		File attributeFile;
		// Creates new or saves over the attribute config file.
		if (fileName == null) {
			attributeFile = new File("AttributesConfig.txt");
		} else {
			attributeFile = new File(fileName);
		}
		attributeFile.createNewFile();
		FileWriter writer = new FileWriter(attributeFile);
		// Short formatting warning
		writer.write("READ: Program requires that each config option has one space, and only one space,\n");
		writer.write("after the equal sign.  If you are editing the options from here, please remember this.\n");
		// Writes all the fields to the text file in the most minimal manner
		// possible.
		for (int i = 0; i < field.length; i++) {
			try {
				writer.write(field[i].getName() + " = " + field[i].get(this)
						+ "\n");
			} catch (Exception e) {
				continue;
			}
		}
		writer.close();
	}
}
