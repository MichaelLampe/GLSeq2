package org.glbrc.glseq2.Project_Files;

public enum Category{
  DATA("Data and Library"),
  REFERENCE("Reference Attributes"),
  RUN("Run Attributes"),
  PREPROCESS("Pre-Processing Attributes"),
  COMMON_PROCESS("Common Processing Options"),
  RSEM("RSEM Attributes"),
  ENVIRONMENT("Environment Attributes");
  
  public final String name;
  Category(String name){
    this.name = name;
  }
}