package org.glbrc.glseq2;

///////////////////////////////////////////////////////////////////////////////
//                   
// Main Class File:  UI.java
// File:             RunOptions.java
//
//
// Author:           Michael Lampe | MrLampe@Wisc.edu
///////////////////////////////////////////////////////////////////////////////

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class Run {
  private String updateFromDatabase;
  private String processedData;
  private String alignment;
  private String counting;
  private String collectResults;
  private String runId;
  private String protocolId;
  private String attributeFilePath;
  private String ampersand;

  /**
   * Instantiates with everything off.
   */
  public Run() {
    updateFromDatabase = "noupdate";
    processedData = "nodataprep";
    alignment = "noalignment";
    counting = "nocounting";
    collectResults = "nocollect";
    runId = "";
    protocolId = "0";
    attributeFilePath = "";
    ampersand = "&";
  }

  /**
   * Converts script into script directory for uses.
   * 
   * @param scriptDirectory
   *          Where the script is located
   * @return the full path directory
   */
  public String getScriptDirectory(String scriptDirectory) {
    File attributeFile = null;
    Scanner scnr = null;
    try {
      attributeFile = new File(attributeFilePath);
      if (attributeFile.equals(null)) {
        scnr = new Scanner(attributeFile);
        while (scnr.hasNextLine()) {
          String[] temp = scnr.nextLine().split(" <- ");
          System.out.println(temp[0]);
          if (temp[0].equals("base.dir")) {
            String returnString = temp[1].replace("\"", "");
            scnr.close();
            return returnString;
          }
        }
      }
    } catch (Exception e) {
      e.printStackTrace(System.out);
      return scriptDirectory;
    }
    return scriptDirectory;
  }

  /**
   * Converts input to usable form.
   * 
   * @param destinationDirectory
   *          Where the script is going
   * @return Converted destination to usable format
   */
  public String getDestinationDirectory(String destinationDirectory) {
    File attributeFile = new File(attributeFilePath);
    Scanner scnr = null;
    try {
      scnr = new Scanner(attributeFile);
      while (scnr.hasNextLine()) {
        String[] temp = scnr.nextLine().split(" <- ");
        if (temp[0].equals("dest.dir.base")) {
          String returnString = temp[1].replace("\"", "");
          scnr.close();
          return returnString;
        }
      }
    } catch (Exception e) {
      System.out.println("Error in loading attribute file");
      System.out.println("Using most recent attribute.");
      return destinationDirectory;
    }
    scnr.close();
    return destinationDirectory;
  }

  /**
   * Gives back a list of all the args for the script to use.
   * 
   * @return List of arguements
   */
  public List<String> returnArgs() {
    List<String> args = new ArrayList<String>();
    args.add("Rscript");
    args.add("GLSeq.top.R");
    args.add(updateFromDatabase);
    args.add(processedData);
    args.add(alignment);
    args.add(counting);
    args.add(collectResults);
    args.add(runId);
    args.add(protocolId);
    args.add(attributeFilePath);
    args.add(ampersand);
    return args;

  }

  // #######################################################
  // ################# Update from Config ##################
  // #######################################################
  /**
   * Updates settings from a configuration file.
   */
  public void updateFromConfig() {
    Field field = null;
    Scanner scnr = null;
    File config = null;
    config = new File("RunConfig.txt");
    try {
      scnr = new Scanner(config);
    } catch (FileNotFoundException e) {
      System.out.println("Run Config file not found");
    }
    while (scnr.hasNextLine()) {
      String temp = scnr.nextLine();
      if (temp.contains("=")) {
        String[] tempArray = temp.split(" = ");
        try {
          field = this.getClass().getDeclaredField(tempArray[0]);
        } catch (NoSuchFieldException e) {
          // Catches if the field does not exist but is in the file passed in
        } catch (SecurityException e) {
          // Fine
        }
        try {
          field.set(this, tempArray[1]);
        } catch (IllegalArgumentException e) {
          // Fine
        } catch (IllegalAccessException e) {
          // Fine
        } catch (ArrayIndexOutOfBoundsException e) {
          // Fine
        }
      }
    }
    scnr.close();
  }

  // #######################################################
  // ################# Saves Current Configs ###############
  // #######################################################
  /** Saves the current fields as a configuration file
   * 
   * @param fileName - What to name it
   * @throws IOException .
   */
  public void saveConfigFile(String fileName) throws IOException {
    final Field[] field = this.getClass().getDeclaredFields();
    File runFile;
    if (fileName == null) {
      runFile = new File("RunConfig.txt");
    } else {
      runFile = new File(fileName);
    }
    runFile.createNewFile();
    FileWriter writer = new FileWriter(runFile);
    writer
        .write("READ: Program requires that each config "
            + "option has one space, and only one space,\n");
    writer
        .write("after the equal sign.  If you are editing"
            + " the options from here, please remember this.\n\n");
    for (int i = 0; i < field.length; i++) {
      try {
        writer.write(field[i].getName() + " = " + field[i].get(this) + "\n");
      } catch (Exception e) {
        // I know this doesn't look as good as it could, but it is just
        // a way to make sure it gets written even if the user has weird
        // stuff in the file. When reading the program won't care about
        // that so it should be fine. (#FamousLastWords)
        continue;
      }
    }
    writer.close();
  }

  // #######################################################
  // ################# Getters ############################
  // #######################################################
  public String getAttributeFilePath() {
    return attributeFilePath;
  }

  public String getProtocolId() {
    return protocolId;
  }

  public String getRunId() {
    return runId;
  }

  public String getCollectResults() {
    return collectResults;
  }

  public String getAlignment() {
    return alignment;
  }

  public String getCounting() {
    return counting;
  }

  public String getProcessedData() {
    return processedData;
  }

  public String getUpdateFromDatabase() {
    return updateFromDatabase;
  }

  public String getAmpersand() {
    return ampersand;
  }

  // #######################################################
  // ################# Setters #############################
  // #######################################################
  public void setAttributeFilePath(String attributeFilePath) {
    this.attributeFilePath = attributeFilePath;
  }

  public void setProtocolId(String protocolId) {
    this.protocolId = protocolId;
  }

  public void setRunId(String runId) {
    runId = runId.replaceAll(" ", "");
    this.runId = runId;
  }

  public void setCollectResults(String collectResults) {
    this.collectResults = collectResults;
  }

  public void setAlignment(String alignment) {
    this.alignment = alignment;
  }

  public void setCounting(String counting) {
    this.counting = counting;
  }

  public void setProcessedData(String processedData) {
    this.processedData = processedData;
  }

  public void setUpdateFromDatabase(String updateFromDatabase) {
    this.updateFromDatabase = updateFromDatabase;
  }

  public void setAmpersand(String ampersand) {
    this.ampersand = ampersand;
  }
}
