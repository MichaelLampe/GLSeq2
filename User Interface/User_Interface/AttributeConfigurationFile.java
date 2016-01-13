package application;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class AttributeConfigurationFile extends AttributeActions {
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

}
