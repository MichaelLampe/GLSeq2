package application;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

public class AttributeAndConfigFileHandler {
	/**
	 * 
	 * @param fileName
	 *            The name of the file the user wishes to create.
	 * @throws IOException .
	 */
	public void saveConfigFile(String fileName) throws IOException {
		/*
		 * An array of all the fields
		 */
		ArrayList<AttributeFactorySingleton.Attribute> attributeArray = AttributeFactorySingleton
				.getInstance().getAllAttributes();

		/*
		 * Create a new attribute config file
		 */
		File attributeFile;
		if (fileName == null) {
			attributeFile = new File("AttributesConfig.txt");
		} else {
			attributeFile = new File(fileName);
		}
		
		/*
		 * Create the file and a file writer to write the upcoming lines.
		 */
		attributeFile.createNewFile();

		// Try w/ resource
		try (FileWriter writer = new FileWriter(attributeFile)) {
			
			/*
			 * Short formatting warning indicating the file format.
			 */
			writer.write("READ: Program requires that each config option has one space, and only one space,\n");
			writer.write("after the equal sign.  If you are editing the options from here, please remember this.\n");
			
			/*
			 * Writes all the fields to the text file in the most minimal manner
			 * possible.
			 */
			for (AttributeFactorySingleton.Attribute attribute : attributeArray) {
				writer.write(attribute.getName() + " = " + attribute.getValue()
						+ "\n");
			}
		}
	}
	
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

		/*
		 * Tracks if attributes were matched to a value or not (If within this
		 * list, it was matched)
		 */
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

				// Variable assignment in R
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
		updateUnmatchedAttributesWithDefault(matchedAttributes);
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

		/*
		 * Tracks if attributes were matched to a value or not (If within this
		 * list, it was matched)
		 */
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
		updateUnmatchedAttributesWithDefault(matchedAttributes);
	}

	/*
	 * Takes care of updating any unmatched fields from the config files with
	 * their default values.
	 */
	private void updateUnmatchedAttributesWithDefault(
			ArrayList<AttributeFactorySingleton.Attribute> matchedAttributes) {
		ArrayList<AttributeFactorySingleton.Attribute> attributeArray = AttributeFactorySingleton
				.getInstance().getAllAttributes();

		/*
		 * For every value that didn't have a matched field, we just set the
		 * value to the default.
		 */
		for (AttributeFactorySingleton.Attribute attribute : attributeArray) {
			if (!matchedAttributes.contains(attribute)) {
				attribute.setValue(attribute.getDefault());
			}
		}
		UpdateUserInterfaceSingleton.getInstance().updateDefaults();
	}
}
