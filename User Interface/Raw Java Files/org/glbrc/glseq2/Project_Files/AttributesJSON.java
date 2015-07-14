package org.glbrc.glseq2.Project_Files;

public enum AttributesJSON {
  directory(Category.DATA.name, "The location containing the raw files."), unzipped(
      Category.DATA.name, "\"TRUE\",\"FALSE\"", "TRUE",
      "If the unprocessed data files are zipped (.gz) or not (.fq or .fastq)"), directoryFq(
      Category.DATA.name, "Directory containing ready-to-go fq or fastq files."), rawFileNames(
      Category.DATA.name, "The raw files names that you would like to process"), strain(
      Category.DATA.name, "-DUMMMY-"), pairedEnd(Category.DATA.name, "\"TRUE\",\"FALSE\"", "TRUE",
      "If the data is paired or single ended"), seqPlatform(Category.DATA.name,
      "\"illumina\",\"capillary\", \"ls454\", \"solid\", \"helicos\", \"iontorrent\",\"pacbio\"",
      "illumina", "The sequencing machine that was used"), qualityScores(Category.DATA.name,
      "\"phred33\",\"phred64\"", "phread33",
      "How the quality of base identification was calculated"), libstrands(Category.DATA.name,
      "\"NULL\",\"F\",\"R\"", "R", "The strandedness of the library"), libNchar(Category.DATA.name,
      "The number of unique characters at the beginning of each file's name"), libList(
      Category.DATA.name, "NULL", "The subset of libraries to process"), countableSamDir(
      Category.DATA.name, "The directory containing already aligned SAM files that can be counted."), presplit(
      Category.DATA.name,
      "\"TRUE\",\"FALSE\"",
      "TRUE",
      "If the paired ended data is split into two files or not.  If your data is single ended you can ignore this option."), prevRunDirectory(
      Category.DATA.name,
      "If just collecting data, this is the previous folder that was created when the run started"), prevRunName(
      Category.DATA.name,
      "If collecting data from a previous run, this is the name that was used for that run."),

  referenceGenome(
      Category.REFERENCE.name,
      "The reference genome will be used to index certain methods such as RSEM or BWA. Should be located in the script folder. The name contains no file extension."), referenceFasta(
      Category.REFERENCE.name, "Name of the reference FASTA file."), referenceGff(
      Category.REFERENCE.name,
      "Name of the reference genomic features file.  Gtf or Gff files are able to be used here. Located in the same folder as the reference genome."), gtfFeatureColumn(
      Category.REFERENCE.name, "9",
      "The number of columns in the GTF file with the gene / other IDs character string"), idAttr(
      Category.REFERENCE.name, "gene_id", "GFF attribute to be used as feature ID"),

  destinationDirectory(
      Category.RUN.name,
      "The directory where a folder will be created that contains all of the results and materials used in this run."), featureCounts(
      Category.RUN.name, "\"Feature Counts\",\"\"", "",
      "Select to run the Feature Counts counting protocol"), rsem(Category.RUN.name,
      "\"RSEM\",\"\"", "",
      "Select to run the RSEM counting protocol. This protocol does not work with gapped aligners."), htseq(
      Category.RUN.name, "\"HTSeq\",\"\"", "", "Select to run the HTSeq counting protocol"), cufflinks(
      Category.RUN.name, "\"Cufflinks\",\"\"", "", "Select to run the Cufflinks counting protocol"), alignmentAlgor(
      Category.RUN.name, "\"BWA\",\"Bowtie\",\"Bowtie2\",\"Cushaw\",\"Cushaw-GPU\",\"TopHat\"",
      "BWA", "The aligner you would like to use to align your reads."), numberCores(
      Category.RUN.name, "The number of cores that your alignment should use."), numberStreams(
      Category.RUN.name, "The number of parallel alignment streams that should happen"), numberStreamsDataPrep(
      Category.RUN.name, "The number of parallel data preparation streams that should happen"), gpuAcceleration(
      Category.RUN.name,
      "\"TRUE\",\"FALSE\"",
      "TRUE",
      "Utilizes a GPU (Available on Scarcity Server 10) to accelerate alignment speed.  Only available when using CUSHAW-GPU"),

  readTrim(Category.PREPROCESS.name, "\"TRUE\",\"FALSE\"", "FALSE",
      "Trim reads during preprocessing. This removes some parts of the sequence prior to alignment."), trimHead(
      Category.PREPROCESS.name, "0", "Automatic headcrop for trimmomatic"), minTrim(
      Category.PREPROCESS.name, "0", "The minimum length of a trimmed read."), artificialFasta(
      Category.PREPROCESS.name,
      "Name of the FASTA file with artificial sequences (Adapters, primers, etc.).  Located in the script directory."),

  strandExtract(Category.COMMON_PROCESS.name, "\"TRUE\",\"FALSE\"", "FALSE",
      "Extracts explicit forward and reverse coverage from the original BAM file."),

  compConf(Category.RSEM.name, "\"TRUE\",\"FALSE\"", "FALSE", "Compute RSEM confidence intervals?"), fragMaxLength(
      Category.RSEM.name, "1000", "Maximal length of a fragment in a paired-end library."), ciMem(
      Category.RSEM.name, "4096",
      "Maximum size (MB) of the auxiliary buffer used to compute confidence intervals"), genobam(
      Category.RSEM.name, "\"TRUE\",\"FALSE\"", "FALSE", "Output a Genome BAM file?"),

  scriptDirectory(Category.ENVIRONMENT.name, "The location of the R Script Files. (GLSeq.Top.R)"), trimPath(
      Category.ENVIRONMENT.name, "Path to teh Trimmomatic"), picardToolsPath(
      Category.ENVIRONMENT.name, "Path to the PIcardTools jar directory"), fastqcPath(
      Category.ENVIRONMENT.name, "Path to the FastQc"), bwaPath(Category.ENVIRONMENT.name,
      "/opt/bifxapps/bin/bwa", "Path to BWA Aligner"), bam2wigPath(Category.ENVIRONMENT.name,
      "Path to the shell script that converts Bam to Wig"), cushawPath(Category.ENVIRONMENT.name,
      "/opt/bifxapps/cushaw2-v2.1.11/cushaw2", "Path to the CUSHAW Aligner"), cushawGpuPath(
      Category.ENVIRONMENT.name, "/opt/bifxapps/cushaw2-gpu-2.1.8-r16/cushaw2-gpu",
      "Path to the Cushaw-GPU Aligner"), cushawIndexPath(Category.ENVIRONMENT.name,
      "/opt/bifxapps/cushaw2-v2.1.11/cushaw2_index/cushaw2_index", "Path to the CUSHAW Indexer"), topHatPath(
      Category.ENVIRONMENT.name, "Path to the Tophat Aligner"), destDirTest(
      Category.ENVIRONMENT.name,
      "Path to a logging file.  This will create general logs and put them in this directory.");

  public final String category;
  public final String options;
  public final String defaultVal;
  public final String description;

  public static final int size = AttributesJSON.values().length;
  
  AttributesJSON(String category, String description) {
    this.category = category;
    this.description = description;
    options = null;
    defaultVal = null;
  }

  AttributesJSON(String category, String defaultVal, String description) {
    this.category = category;
    this.description = description;
    this.defaultVal = defaultVal;
    options = null;
  }

  AttributesJSON(String category, String options, String defaultVal, String description) {
    this.category = category;
    this.options = options;
    this.defaultVal = defaultVal;
    this.description = description;
  }
}