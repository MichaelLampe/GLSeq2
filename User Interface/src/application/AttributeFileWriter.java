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
import java.util.Date;

public class AttributeFileWriter {
	/**
	 * 
	 * @return AttributeFileLocation - The new attribute file's path
	 * @throws IOException
	 *             .
	 */
	public String writeAttributesFile(String file_name, String file_location) throws IOException {
		/*
		 * Writes a new attribute file based on the current format that we have
		 * outlined. Saves it as an R file, with the current date. The file to
		 * be created
		 */
		File attributeFile = null;

		// Adds the date to the file
		Date current = new Date();
		SimpleDateFormat day = new SimpleDateFormat("yyyy.MM.dd");
		SimpleDateFormat time = new SimpleDateFormat("HH.mm.ss");
		String date_time = "_" + day.format(current) + "_" + time.format(current);

		/*
		 * If null, names file for user
		 */
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
		String[] lines = fullFileAsString().split("\n");
		for (String line : lines) {
			line = line + "\n";
			writer.write(line);
		}
		writer.close();
		return attFileLocation;
	}

	public String fullFileAsString() {
		StringBuilder build = new StringBuilder("");

		String header = "# Great Lakes Seq package for low-level processing of RNA-Seq data\n # \n";
		build.append(header);

		for (AttributeFactorySingleton.Attribute attribute : AttributeFactorySingleton.getInstance()
				.getAllAttributes()) {
			/*
			 * Script attributes are not used in attribute file, but in other
			 * instances
			 */
			if (!attribute.getCategory().equals("script")) {
				build.append(attributeFileGenerateLine(attribute));
			}
		}
		String footer = "# Collect all the counting algos together \n";
		build.append(footer);

		build.append("cAlgor <- c(HTSeq,RSEM,FeatureCounts,Cufflinks,rockhopperCount)\n");

		/*
		 * @ DEPRECIATED: We depreciated these two, but it is easier to make
		 * sure they still run if the user has an older version of the Rscript.
		 * Very unlikely, but doesn't hurt.
		 */
		build.append("raw.dir <- data.directory\n");
		build.append("readyData.dir <- data.directory\n");

		return build.toString();
	}

	private String attributeFileGenerateLine(AttributeFactorySingleton.Attribute attribute) {
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
