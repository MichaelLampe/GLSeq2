package org.glbrc.glseq2;

public final class ButtonEnums {
  public enum AttributeButton {
    UPDATE("Updating from Database"), NO_UPDATE("NOT Updating from Database"),
    //
    PREPROCESSING("Pre-Processing Data"), NO_PREPROCESSING("NOT Pre-Processing Data"),
    //
    ALIGNMENT("Aligning"), NO_ALIGNMENT("NOT Aligning"),
    //
    COUNT("Counting"), NO_COUNT("NOT Counting"),
    //
    COLLECT("Collecting Results"), NO_COLLECT("NOT Collecting Results"),
    //
    BACKGROUND("Running in the Background"), NO_BACKGROUND("NOT Running in the Background");

    public final String value;

    private AttributeButton(String value) {
      this.value = value;
    }
  }

  public enum OptionButton {
    // Data and Library
    ZIPPED("Using Zipped Files (.gz)"), UNZIPPED("Using Unzipped Files (.fq)"), SINGLE(
        "Using Single Ended Data"), PAIRED("Using Paired Ended Data"),
    // Processing
    TRIMMING("Trimming Raw Reads"), NO_TRIMMING("Not Trimming Raw Reads"),
    //
    // Script Running Options
    EXTRACT("Extracting Forward and Reverse Coverage from Original BAM File"), 
    NO_EXTRACT("NOT Extracting Forward and Reverse Coverage from Original BAM File"), 
    //
    COMPUTE("Computing Confidence Intervals"), 
    NO_COMPUTE("NOT Computing Confidence Intervals"), 
    //
    OUTPUT("Outputting Genome BAM"),
    NO_OUTPUT("NOT Outputting Genome Bam"), 
    //
    PRESPLIT("Presplit"),
    NO_PRESPLIT("NOT Presplit");

    public final String value;

    private OptionButton(String value) {
      this.value = value;
    }
  }

  public enum Attribute {
    /*
     * These are all based on whatever the arguments that are requested by the
     * Rscript are
     */
    UPDATE("update"), NO_UPDATE("noupdate"),
    //
    PREPROCESSING("dataprep"), NO_PREPROCESSING("nodataprep"),
    //
    ALIGNMENT("alignment"), NO_ALIGNMENT("noalignment"),
    //
    COUNT("counting"), NO_COUNT("nocounting"),
    //
    COLLECT("collect"), NO_COLLECT("nocollect"),
    //
    BATCH("Running all queued runs"), NO_BATCH("Run only selected queued run");

    public final String value;

    private Attribute(String value) {
      this.value = value;
    }
  }
}
