package application;

// This contains the name (Enum name), default values, possible values, and tool tips.
public enum AttributesJSON {
  dataDirectory("data.directory",Category.DATA.name, "The location containing the data files."),
  
  storageDestination("storage.destination",Category.DATA.name, "The place where the files will be finally stored"),
  
  unzipped("unzipped",Category.DATA.name, "\"TRUE\",\"FALSE\"", "TRUE",
      "If the unprocessed data files are zipped (.gz) or not (.fq or .fastq)"), 
  strain("strain",
      Category.DATA.name, "Strain of data"), 
  pairedEnd("paired.end",Category.DATA.name, "\"TRUE\",\"FALSE\"", "TRUE",
      "If the data is paired or single ended"), 
  seqPlatform("seqPlatform",Category.DATA.name,
      "\"illumina\",\"capillary\", \"ls454\", \"solid\", \"helicos\", \"iontorrent\",\"pacbio\"",
      "illumina", "The sequencing machine that was used"), 
  qualityScores("qScores",Category.DATA.name,
      "\"phred33\",\"phred64\"", "phread33",
      "How the quality of base identification was calculated"), 
  libStrands("libstrand",Category.DATA.name,
      "\"NULL\",\"F\",\"R\"", "R", "The strandedness of the library"), 
  libNchar("libNchar",Category.DATA.name,
      "The number of unique characters at the beginning of each file's name"), 
  countableSamDir("countable.sams.dir",
      Category.DATA.name, "The directory containing already aligned SAM files that can be counted."), 
  presplit("presplit",
      Category.DATA.name,
      "\"TRUE\",\"FALSE\"",
      "TRUE",
      "If the paired ended data is split into two files or not.  If your data is single ended you can ignore this option."), 
  prevRunDirectory("previous.dir",
      Category.DATA.name,
      "If just collecting data, this is the previous folder that was created when the run started"), 
  prevRunName("previous.run.name",
      Category.DATA.name,
      "If collecting data from a previous run, this is the name that was used for that run."),
  referenceGenome("rGenome",
      Category.REFERENCE.name,
      "The reference genome will be used to index certain methods such as RSEM or BWA. Should be located in the script folder. The name contains no file extension."), 
  referenceFasta("refFASTAname",
      Category.REFERENCE.name, "Name of the reference FASTA file."), 
  referenceGff("refGFFname",
      Category.REFERENCE.name,
      "Name of the reference genomic features file.  Gtf or Gff files are able to be used here. Located in the same folder as the reference genome."), 
  gtfFeatureColumn("gtfFeatureColumn",
      Category.REFERENCE.name, "9",
      "The number of columns in the GTF file with the gene / other IDs character string"), 
  idAttr("idAttr",
      Category.REFERENCE.name, "gene_id", "GFF attribute to be used as feature ID"),

  destinationDirectory("dest.dir.base",
      Category.RUN.name,
      "The directory where a folder will be created that contains all of the results and materials used in this run."), 
  FeatureCounts("FeatureCounts",
      Category.RUN.name, "\"Feature Counts\",\"\"", "",
      "Select to run the Feature Counts counting protocol"), 
  RSEM("RSEM",Category.RUN.name,
      "\"RSEM\",\"\"", "",
      "Select to run the RSEM counting protocol. This protocol does not work with gapped aligners."), 
  HTseq("HTSeq",
      Category.RUN.name, "\"HTseq\",\"\"", "", "Select to run the HTSeq counting protocol"), 
      
      
  rockhopperCount("Rockhopper",
      Category.RUN.name, "\"Rockhopper\",\"\"", "", "Select to run the Rockhopper counting protocol"), 
      
      
  Cufflinks("Cufflinks",
      Category.RUN.name, "\"Cufflinks\",\"\"", "", "Select to run the Cufflinks counting protocol"), 
  aAlgor("aAlgor",
      Category.RUN.name, "\"BWA\",\"Bowtie\",\"Bowtie2\",\"Cushaw\",\"Cushaw_GPU\",\"TopHat\",\"Rockhopper\"",
      "BWA", "The aligner you would like to use to align your reads."), 
  numberCores("nCores",
      Category.RUN.name, "The number of cores that your alignment should use."), 
  numberStreams("nStreams",
      Category.RUN.name, "The number of parallel alignment streams that should happen"), 
  numberStreamsDataPrep("nStreamsDataPrep",
      Category.RUN.name, "The number of parallel data preparation streams that should happen"),

  readTrim("readTrim",Category.PREPROCESS.name, "\"TRUE\",\"FALSE\"", "FALSE",
      "  reads during preprocessing. This removes some parts of the sequence prior to alignment."), 
  trimHead("trimhead",
      Category.PREPROCESS.name, "0", "Automatic headcrop for trimmomatic"), 
  minTrim("minTrim",
      Category.PREPROCESS.name, "0", "The minimum length of a trimmed read."), 
  artificialFasta("artificial.fq",
      Category.PREPROCESS.name,
      "Name of the FASTA file with artificial sequences (Adapters, primers, etc.).  Located in the script directory."),

  fragMaxLength("fragMaxLength",
      Category.RSEM.name, "1000", "Maximal length of a fragment in a paired-end library."),
  ciMem("ciMem",
      Category.RSEM.name, "4096",
      "Maximum size (MB) of the auxiliary buffer used to compute confidence intervals"), 

  scriptDirectory("base.dir",Category.ENVIRONMENT.name, "The location of the R Script Files. (GLSeq.Top.R)"), 
  trimPath("trimPath",
      Category.ENVIRONMENT.name,"/opt/bifxapps/Trimmomatic-0.30/trimmomatic-0.30.jar", "Path to the Trimmomatic"), 
  picardToolsPath("picardToolsPath",
      Category.ENVIRONMENT.name,"/opt/bifxapps/picard-tools/","Path to the PIcardTools jar directory"), 
  fastqcPath("fastqcPath",
      Category.ENVIRONMENT.name,"/opt/bifxapps/bin/fastqc","Path to the FastQc"), 
  bwaPath("bwaPath",Category.ENVIRONMENT.name,
      "/opt/bifxapps/bin/bwa", "Path to BWA Aligner"), 
  bam2wigPath("bam2wigPath",
      Category.ENVIRONMENT.name,"/home/GLBRCORG/omoskvin/run/bam2wig.sh",
      "Path to the shell script that converts Bam to Wig"), 
  cushawPath("CUSHAW.path",
      Category.ENVIRONMENT.name,
      "/opt/bifxapps/cushaw2-v2.1.11/cushaw2", "Path to the CUSHAW Aligner"), 
  cushawGpuPath("CUSHAW.GPU.path",
      Category.ENVIRONMENT.name, "/opt/bifxapps/cushaw2-gpu-2.1.8-r16/cushaw2-gpu",
      "Path to the Cushaw-GPU Aligner"), 
  cushawIndexPath("CUSHAW.index.path",
      Category.ENVIRONMENT.name,
      "/opt/bifxapps/cushaw2-v2.1.11/cushaw2_index/cushaw2_index", "Path to the CUSHAW Indexer"),
  topHatPath("TopHat.path",
      Category.ENVIRONMENT.name,"/opt/bifxapps/tophat-2.0.11.Linux_x86_64/tophat2", "Path to the Tophat Aligner"),
  rockhopperPath("Rockhopper.path",
          Category.ENVIRONMENT.name, "Path to the Rockhopper Aligner"),
  destDirTest("destDirTest",
      Category.ENVIRONMENT.name,
      "Path to a logging file.  This will create general logs and put them in this directory."),
  starPath("star.path",
      Category.ENVIRONMENT.name,
      "Path to a logging file.  This will create general logs and put them in this directory."),
  hisatPath("hisat.path",Category.ENVIRONMENT.name,
      "Path to the folder containing the hisat-build and hisat files."),    
      
  uiDivider("uiDivider",
          Category.SCRIPT.name,
          "Divider for Glow"),
  Update("Updating",Category.SCRIPT.name,"\"noupdate\",\"update\"","noupdate",
      "If to update from glow."),
  DataPrep("Data Preparation",Category.SCRIPT.name,"\"dataprep\",\"nodataprep\"","nodataprep",
      "If to prepare raw data."),
  Alignment("Alignment",Category.SCRIPT.name,"\"alignment\",\"noalignment\"", "noalignment",
      "If to align data."),
  Counting("Counting",Category.SCRIPT.name,"\"counting\",\"nocounting\"", "nocounting",
      "If to count data."),
  Collect("Collect",
      Category.SCRIPT.name,"\"collect\",\"nocollect\"", "nocollect",
      "If to collect results summary."),
  RunName("Run Name",
      Category.SCRIPT.name,
      "The name of the run to be started.");

  public final String category;
  public final String options;
  public final String defaultVal;
  public final String description;
  public final String field_name;

  public static final int size = AttributesJSON.values().length;
  
  AttributesJSON(String name,String category, String description) {
    this.field_name = name;
    this.category = category;
    this.description = description;
    options = null;
    defaultVal = "NULL";
  }

  AttributesJSON(String name,String category, String defaultVal, String description) {
    this.field_name = name;
    this.category = category;
    this.description = description;
    this.defaultVal = defaultVal;
    options = null;
  }

  AttributesJSON(String name,String category, String options, String defaultVal, String description) {
    this.field_name = name;
    this.category = category;
    this.options = options;
    this.defaultVal = defaultVal;
    this.description = description;
  }
}