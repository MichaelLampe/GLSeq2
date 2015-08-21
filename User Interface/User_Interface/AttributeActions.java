package application;

///////////////////////////////////////////////////////////////////////////////
//                   
// Main Class File:  UI.java
// File:             Attributes.java  
//
//
// Author:           Michael Lampe | MrLampe@Wisc.edu
///////////////////////////////////////////////////////////////////////////////
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Scanner;

public final class AttributeActions {

  /**
   * Creates a new or replaces the current attribute file and adds all the field
   * data to that attribute file.
   * 
   *
   */
  public void setAttributes() {
    Scanner scnr = null;
    File config = null;
    try {
      config = new File("AttributesConfig.txt");
      scnr = new Scanner(config);
    } catch (Exception e) {
      return;
      // This is expected
      // If there is no AttributeFile because of no previous runs
    }
    // Looks through the whole document for variables that match field names
    while (scnr.hasNextLine()) {
      String temp = scnr.nextLine();
      if (temp.contains("=")) {
        String[] tempArray = temp.split(" = ");
        try {
          Attributes.getInstance().attributesCollection.get(tempArray[0]).setValue(tempArray[1]);
        } catch (Exception e) {
          // No matching value
        }
      }
    }
    scnr.close();
    UpdateUI.getInstance().updateDefaults();
  }

  /**
   * Creates a new or replaces the current attribute file and adds all the field
   * data to that attribute file.
   * 
   * <p>
   * This method uses an input attributeFile.R to take old data and load it in.
   *
   */
  public void setAttributes(File attributeFile) {
    ArrayList<Attribute> matchedAttributes = new ArrayList<Attribute>();
    Scanner scnr = null;
    try {
      scnr = new Scanner(attributeFile);
      // Looks through the whole document for variables that match field
      // names
      while (scnr.hasNextLine()) {
        String temp = scnr.nextLine();
        if (temp.contains("<-")) {
          String[] tempArray = temp.split(" <- ");
          try {
            tempArray[1] = tempArray[1].replace("\"","");
            tempArray[1] = tempArray[1].replace("\'","");
            Attributes.getInstance().attributesCollection.get(tempArray[0]).setValue(tempArray[1]);
            matchedAttributes.add(Attributes.getInstance().attributesCollection.get(tempArray[0]));
          } catch (Exception e) {
            // if doesn't map
          }
        }
      }
    } catch (Exception e) {
      System.out.println("Attribute file not found");
    } finally {
      ArrayList<Attribute> attributeArray = new ArrayList<Attribute>(
          Attributes.getInstance().attributesCollection.values());
      for (Attribute attribute : attributeArray) {
        if (!matchedAttributes.contains(attribute)){
          attribute.setValue(attribute.getDefault());
        }
      }
      // Makes sure to close the scanner
      scnr.close();
      UpdateUI.getInstance().updateDefaults();
    }
  }

  /**
   * Future use case for possible SSH construction of AttributeFiles in a
   * programmatic manner. Otherwise also let's the user not use the UI when
   * accessing it from a server.
   * 
   * @param inputArgs
   *          An array of strings passed in from another source
   */
  public void setAttributes(String[] inputArgs) {
    ArrayList<Attribute> attributeArray = new ArrayList<Attribute>(
        Attributes.getInstance().attributesCollection.values());
    for (String param : inputArgs) {
      String[] parts = param.split("=");
      // Complexity is fine when it is cheap.
      for (Attribute att : attributeArray) {
        if (att.getUiName().equals(parts[0])) {
          Attribute current = Attributes.getInstance().attributesCollection.get(att.getName());
          try {
            current.setValue(parts[1]);
          } catch (Exception e) {
            // Improper field was entered
          }
        }
      }
    }
  }

  /**
   * 
   * @param fileName
   *          The name of the file the user wishes to create.
   * @throws IOException .
   */
  public void saveConfigFile(String fileName) throws IOException {
    // An array of all the fields
    ArrayList<Attribute> attributeArray = new ArrayList<Attribute>(
        Attributes.getInstance().attributesCollection.values());
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
    writer.write("READ: Program requires that each config "
        + "option has one space, and only one space,\n");
    writer.write("after the equal sign.  If you are editing "
        + "the options from here, please remember this.\n");
    // Writes all the fields to the text file in the most minimal manner
    // possible.
    for (Attribute attribute : attributeArray) {
      writer.write(attribute.getName() + " = " + attribute.getValue() + "\n");
    }
    writer.close();
  }

  private String checkSlash(String file_name) {
    if (file_name.substring(file_name.length() - 1).equals("/")) {
      return file_name;
    } else {
      file_name = file_name + "/";
      return file_name;
    }
  }

  public void returnJson() {
    System.out.println("{\"attributes\":[");
    ArrayList<Attribute> attributeArray = new ArrayList<Attribute>(
        Attributes.getInstance().attributesCollection.values());
    for (int i = attributeArray.size() - 1; i >= 0; i--) {
      if (attributeArray.get(i).getCategory().equals(Category.SCRIPT.name)) {
        attributeArray.remove(i);
      }
    }
    for (int i = 0; i < attributeArray.size(); i++) {
      System.out.print(formatJson(attributeArray.get(i).getUiName()));
      if (i < attributeArray.size() - 1) {
        System.out.print(",");
      }
      System.out.println("");
    }
    System.out.println("]}");
  }

  /*
   * JSON formatted String to return to GLOW
   */
  private String formatJson(String field_name) {
    String jsonKey = "{\"key\":";
    AttributesJSON current = Enum.valueOf(AttributesJSON.class, field_name);
    jsonKey += "\"" + current.name() + "\"";
    if (current.category != null) {
      jsonKey += ", \"category\":\"" + current.category + "\"";
    }
    if (current.options != null) {
      jsonKey += ", \"options\":[" + current.options + "]";
    }
    if (current.defaultVal != null) {
      jsonKey += ", \"default\":\"" + defaultNullIsEmpty(current.defaultVal) + "\"";
    }
    if (current.description != null) {
      jsonKey += ", \"description\":\"" + current.description + "\"";
    }
    jsonKey += "}";
    return jsonKey;
  }

  private String defaultNullIsEmpty(String defaultVal) {
    if (defaultVal.equals("NULL")) {
      return "";
    }
    return defaultVal;
  }

  /**
   * 
   * @return AttributeFileLocation - The new attribute file's path
   * @throws IOException .
   */
  public String writeAttributesFile(String file_name, String file_location) throws IOException {
    // Writes a new attribute file based on the current format that we have
    // outlined.
    // Saves it as an R file, with the current date.
    // The file to be created
    File attributeFile = null;
    // Adds the date to the file
    Date current = new Date();
    SimpleDateFormat day = new SimpleDateFormat("yyyy.MM.dd");
    SimpleDateFormat time = new SimpleDateFormat("HH.mm.ss");
    String date_time = "_" + day.format(current) + "_" + time.format(current);
    // If null, names file for user
    if (file_name == null) {
      if (file_location == null) {
        attributeFile = new File("Attribute_" + date_time + ".R");
      } else {
        file_location = checkSlash(file_location);
        attributeFile = new File(file_location + "Attribute_" + date_time + ".R");
      }
    } else {
      // Custom file name set by the user.
      if (file_location == null) {
        attributeFile = new File(file_name + date_time + ".R");
      } else {
        file_location = checkSlash(file_location);
        attributeFile = new File(file_location + file_name + date_time + ".R");
      }
    }
    attributeFile.createNewFile();
    final String attFileLocation = attributeFile.getAbsolutePath();
    FileWriter writer = new FileWriter(attributeFile);
    // The entire attribute file in this format
    String header = "# Great Lakes Seq package for low-level processing of RNA-Seq data\n # \n";
    writer.write(header);
    // Converts the hashmap to an arraylist for for eachin'
    for (Attribute attribute : new ArrayList<Attribute>(
        Attributes.getInstance().attributesCollection.values())) {
      if (!attribute.getCategory().equals(Category.SCRIPT.name)) {
        writer.write(generateLine(attribute));
      }
    }
    // Special code to make the R script's life easier.
    // Add anything to further process the data here.
    String footer = "# Collect all the counting algos together \n";
    footer += "cAlgor <- c(HTSeq,RSEM,FeatureCounts,Cufflinks,rockhopperCount)\n";
    footer += "raw.dir <- data.directory\n";
    footer += "readyData.dir <- data.directory\n";
    writer.write(footer);
    writer.close();
    return attFileLocation;
  }

  private String generateLine(Attribute attribute) {
    String line = "";
    // # is a comment in R, so we add the tooltip as a comment
    line += "#" + attribute.getToolTip() + "\n";
    String value = attribute.getValue();
    try {
      Integer.valueOf(value);
    } catch (NumberFormatException e) {
      value = checkValue(value);
    }
    line += attribute.getName() + " <- " + value + "\n";
    return line;
  }

  private String checkValue(String value) {
    switch (value) {
      case "TRUE":
      return value;
      case "FALSE":
      return value;
      case "NULL":
      return value;
      default:
      return ("\"" + value + "\"");
    }
  }
}
