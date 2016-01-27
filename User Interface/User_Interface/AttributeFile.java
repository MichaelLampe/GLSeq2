package application;

import java.io.IOException;

/*
 * Don't refactor the name, as people may utilize it already with GLOW.
 */
public class AttributeFile {

	private static final AttributeActionsCommandLine att = new AttributeActionsCommandLine();

	// Can be launched without the other parts of the program by using
	// -cp as the java option and referring to this as application/AttributeFile
	public static void main(String[] args) {
		generateAttributeFile(args);
		System.exit(0);
	}

	/*
	 * Create a new attribute file based on the args. FILE_NAME= changes the
	 * attribute file name FILE_LOCATION= changes where the file is saved All
	 * other args are field names, = , and then the value that should be placed
	 * in that field. Note: There is currently very little validation that goes
	 * on here, as all the params are strings. May look into this in the future.
	 */
	public static String generateAttributeFile(String[] args) {
		try {
			String file_name = null;
			String file_location = null;
			att.setAttributes(args);
			try {
				if (args[0].contains("FILE_NAME")) {
					file_name = args[0].split("=")[1];
					if (args[1].contains("FILE_LOCATION")) {
						file_location = args[1].split("=")[1];
					}
				} else {
					if (args[0].contains("FILE_LOCATION")) {
						file_location = args[0].split("=")[1];
					}
				}
			} catch (IndexOutOfBoundsException e) {
				// Do nothing here.
			}
			System.out.println("ATTIRUBTE_FILE_PATH="
					+ att.writeAttributesFile(file_name, file_location));

			return file_location;
		} catch (IOException e) {
			System.out
					.println("Error constructing attribute file. Likely unable to create attribute file where the JAR file exists.");
		}
		return null;
	}
}
