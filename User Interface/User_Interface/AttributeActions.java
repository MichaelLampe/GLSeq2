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
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Scanner;

public class AttributeActions {
	/**
	 * Creates a new or replaces the current attribute file and adds all the
	 * field data to that attribute file.
	 * 
	 * This method uses an input attributeFile.R to take old data and load it
	 * in.
	 * 
	 * @throws FileNotFoundException
	 *
	 */
	public void setAttributesAttributeFile(File attributeFile)
			throws FileNotFoundException {
		if (!attributeFile.exists()) {
			throw new FileNotFoundException();
		}

		// Tracks if attributes were matched to a value or not (If within this
		// list, it was matched)
		ArrayList<AttributeFactorySingleton.Attribute> matchedAttributes = new ArrayList<AttributeFactorySingleton.Attribute>();

		try (Scanner scnr = new Scanner(attributeFile)) {
			// Looks through the whole document for variables that match field
			// names
			while (scnr.hasNextLine()) {
				String singleLine = scnr.nextLine();
				// This is the R Attribute file variable assignment expression.

				// Means this line is a comment.
				if (singleLine.startsWith("#")) {
					continue;
				}
				if (singleLine.contains("<-")) {
					String[] tempArray = singleLine.split(" <- ");
					if (tempArray.length > 1) {
						String attributeName = tempArray[0];

						String attributeValue = tempArray[1];

						try {
							attributeValue = attributeValue.replace("\"", "");
							AttributeFactorySingleton.Attribute currentAttribute = AttributeFactorySingleton
									.getInstance().getAttribute(attributeName);
							if (currentAttribute != null) {
								currentAttribute.setValue(attributeValue);
								matchedAttributes.add(currentAttribute);
							}
						} catch (Exception e) {
							continue;
						}
					}
				}
			}
		}
		ArrayList<AttributeFactorySingleton.Attribute> attributeArray = AttributeFactorySingleton
				.getInstance().getAllAttributes();

		// For every value that didn't have a matched field, we just set the
		// value to the default.
		for (AttributeFactorySingleton.Attribute attribute : attributeArray) {
			if (!matchedAttributes.contains(attribute)) {
				attribute.setValue(attribute.getDefault());
			}
		}
		// Update the UI here.
		UpdateUserInterfaceSingleton.getInstance().updateDefaults();
	}

	/**
	 * Creates a new or replaces the current attribute file and adds all the
	 * field data to that attribute file.
	 * 
	 * This method uses a saved text config file.
	 * 
	 * @throws FileNotFoundException
	 *
	 */
	public void setAttributesTextConfig(File attributeFile)
			throws FileNotFoundException {
		if (!attributeFile.exists()) {
			throw new FileNotFoundException();
		}

		// Tracks if attributes were matched to a value or not (If within this
		// list, it was matched)
		ArrayList<AttributeFactorySingleton.Attribute> matchedAttributes = new ArrayList<AttributeFactorySingleton.Attribute>();

		try (Scanner scnr = new Scanner(attributeFile)) {
			// Looks through the whole document for variables that match field
			// names
			while (scnr.hasNextLine()) {
				String singleLine = scnr.nextLine();
				// This is the R Attribute file variable assignment expression.
				if (singleLine.contains("=")) {
					String[] tempArray = singleLine.split(" = ");
					// If the split only has one item, it doesn't have anything
					// worth assigning.
					if (tempArray.length > 1) {
						String attributeName = tempArray[0];
						String attributeValue = tempArray[1];

						try {
							tempArray[1] = tempArray[1].replace("\"", "");
							AttributeFactorySingleton.getInstance()
									.setAttributeValue(attributeName,
											attributeValue);

							matchedAttributes.add(AttributeFactorySingleton
									.getInstance().getAttribute(attributeName));
						} catch (Exception e) {
							// if doesn't map
						}
					}
				}
			}
		}
		ArrayList<AttributeFactorySingleton.Attribute> attributeArray = AttributeFactorySingleton
				.getInstance().getAllAttributes();

		// For every value that didn't have a matched field, we just set the
		// value to the default.
		for (AttributeFactorySingleton.Attribute attribute : attributeArray) {
			if (!matchedAttributes.contains(attribute)) {
				attribute.setValue(attribute.getDefault());
			}
		}
		UpdateUserInterfaceSingleton.getInstance().updateDefaults();
	}

	/**
	 * 
	 * @return AttributeFileLocation - The new attribute file's path
	 * @throws IOException .
	 */
	public String writeAttributesFile(String file_name, String file_location)
			throws IOException {
		// Writes a new attribute file based on the current format that we have
		// outlined.
		// Saves it as an R file, with the current date.
		// The file to be created
		File attributeFile = null;
		// Adds the date to the file
		Date current = new Date();
		SimpleDateFormat day = new SimpleDateFormat("yyyy.MM.dd");
		SimpleDateFormat time = new SimpleDateFormat("HH.mm.ss");
		String date_time = "_" + day.format(current) + "_"
				+ time.format(current);
		// If null, names file for user
		if (file_name == null) {
			if (file_location == null) {
				attributeFile = new File("Attribute_" + date_time + ".R");
			} else {
				file_location = checkSlash(file_location);
				attributeFile = new File(file_location + "Attribute_"
						+ date_time + ".R");
			}
		} else {
			// Custom file name set by the user.
			if (file_location == null) {
				attributeFile = new File(file_name + ".R");
			} else {
				file_location = checkSlash(file_location);
				attributeFile = new File(file_location + file_name + ".R");
			}
		}
		attributeFile.createNewFile();
		final String attFileLocation = attributeFile.getAbsolutePath();
		FileWriter writer = new FileWriter(attributeFile);
		// The entire attribute file in this format
		String header = "# Great Lakes Seq package for low-level processing of RNA-Seq data\n # \n";
		writer.write(header);
		// Converts the hashmap to an arraylist for for eachin'
		for (AttributeFactorySingleton.Attribute attribute : AttributeFactorySingleton
				.getInstance().getAllAttributes()) {
			if (!attribute.getCategory().equals("script")) {
				writer.write(generateLine(attribute));
			}
		}
		/*
		 * Special code to make the R script's life easier. Add anything to
		 * further process the data here.
		 */
		String footer = "# Collect all the counting algos together \n";
		footer += "cAlgor <- c(HTSeq,RSEM,FeatureCounts,Cufflinks,rockhopperCount)\n";

		/*
		 * @ DEPRECIATED: We depreciated these two, but it is easier to make
		 * sure they still run if the user has an older version of the Rscript.
		 * Very unlikely, but doesn't hurt.
		 */
		footer += "raw.dir <- data.directory\n";
		footer += "readyData.dir <- data.directory\n";
		writer.write(footer);
		writer.close();
		return attFileLocation;
	}

	private String generateLine(AttributeFactorySingleton.Attribute attribute) {
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

	private String checkSlash(String file_name) {
		if (file_name.substring(file_name.length() - 1).equals("/")) {
			return file_name;
		} else {
			file_name = file_name + "/";
			return file_name;
		}
	}

	private String checkValue(String value) {
		// Check for null pointers
		if (value == null) {
			return "\"\"";
		}
		// if list don't add the quotes
		if (value.startsWith("c(")) {
			return value;
		}
		switch (value) {
		case "TRUE":
			return value;
		case "FALSE":
			return value;
		case "NULL":
			return value;
		default:
			// Escape a trailing \
			if (value.endsWith("\\")) {
				return ("\"" + value + "\\\"");
			}
			return ("\"" + value + "\"");
		}
	}
}
