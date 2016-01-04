package application;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Map;

public class CombinatorialMutator {

	private static ArrayList<String> mutateFields = new ArrayList<String>();
	private static ArrayList<String> attributeFiles = new ArrayList<String>();

	private static ArrayList<ArrayList<String>> combinations = new ArrayList<ArrayList<String>>();

	public static void main(String[] args) {
		if (args.length < 1) {
			printHelp();
			System.exit(0);
		}

		// Map of all the attributes
		Map<String, Attribute> map = Attributes.getInstance().attributesCollection;

		// Go through all the args and parse them
		for (String arg : args) {

			// Loads an attribute file
			if (arg.contains("attributefile=")) {
				loadAttributeFile(arg.split("attributefile=")[1]);

				// Equals mean that you want to assign a value to something.
			} else if (arg.contains("=")) {
				String[] split = arg.split("=");
				String field_name = split[0];
				String field_value = split[1];

				// Assign value to a given field.
				if (map.containsKey(field_name)) {
					Attribute current_attribute = map.get(field_name);
					current_attribute.setValue(field_value);
				}

				// The -- with the field name is a mutator. If there are any
				// options, they will all be generated in individual attribute
				// files.
			} else if (arg.startsWith("--")) {
				String mutateArg = arg.split("--")[1];
				mutateFields.add(mutateArg);
			}
		}
		createCombinations();
	}

	private static void createCombinations() {

		for (int i = 0; i < mutateFields.size(); i++) {

			String[] options = Attributes.getInstance().attributesCollection.get(mutateFields.get(i)).getOptions();

			if (options != null) {

				ArrayList<String> mutatedVariables = new ArrayList<String>();

				for (int j = 0; j < options.length; j++) {
					// Remove white spaces and add it in the form that
					// AttributeFile.java is used to.
					mutatedVariables.add(mutateFields.get(i) + "=" + options[j].trim());
				}

				// Add all the possible combinations of that variable to the
				// mutation array.
				combinations.add(mutatedVariables);
			}
		}

		createNewAttributeFiles(combinations);

	}

	private static void createNewAttributeFiles(ArrayList<ArrayList<String>> mutationArray) {
		for (int i = 0; i < mutationArray.get(0).size(); i++) {
			String currentArgs = mutationArray.get(0).get(i);
			recursiveCreateAttributeFiles(mutationArray, 1, currentArgs);
		}
	}

	private static void recursiveCreateAttributeFiles(ArrayList<ArrayList<String>> mutationArray, int elementIndex,
			String currentArgs) {
		if (elementIndex == mutationArray.size()) {
			String fileName = currentArgs.replaceAll("\\s+", "_");
			fileName = fileName.replace("=", "_");
			fileName = fileName.replace("\"", "");
			currentArgs = "FILE_NAME=" + fileName + " " + currentArgs.trim();
			// Create file once the recursion is at the end
			String fileLocation = AttributeFile.generateAttributeFile(currentArgs.split(" "));
			attributeFiles.add(fileLocation);
			return;
		}

		// Create copy to avoid bad repeating
		String inputArg = currentArgs;

		for (int i = 0; i < mutationArray.get(elementIndex).size(); i++) {
			currentArgs = inputArg + " " + mutationArray.get(elementIndex).get(i);
			recursiveCreateAttributeFiles(mutationArray, elementIndex + 1, currentArgs);
		}
	}

	/*
	 * Prints a help message when desired.
	 */
	private static void printHelp() {
		System.out
				.println("Call this individual class to generate varying combinations of a given set of conditions.\n");
		System.out.println("Anything after a colon is a space separated command-line argument you can pass.\n\n");
		System.out.println("To load data from a prior attribute file: attributefile=[Absolute_Path_Of_Attribute_File]");
		System.out.println("To individually set an attribute: [Attribute_Name (See list below)]=[AttributeValue]");
		System.out.println(
				"To designate that attribute files should be generated for all combinations of a given value: --[Attribute_Name (See list below)]");
		System.out.println("\n\n\nList of currently supported attributes\n");
		Map<String, Attribute> map = Attributes.getInstance().attributesCollection;
		for (Attribute value : map.values()) {
			System.out.println(value.getName() + " - " + value.getToolTip());
		}
		System.out.println("\n\n\nList of all currently supported attributes with multiple possible values");
		for (Attribute value : map.values()) {
			if (value.getOptions() != null) {
				System.out.println(value.getName() + " - " + value.getToolTip());
				System.out.print("\tPossible values are: ");
				for (String opt : value.getOptions()) {
					System.out.print(opt + " ");
				}
				System.out.println();
			}
		}
	}

	/*
	 * Sets a given value to attribute files.
	 */
	private static boolean loadAttributeFile(String file_name) {
		AttributeActions action = new AttributeActions();
		File att_file = new File(file_name);

		try {
			action.setAttributesAttributeFile(att_file);
			return true;
		} catch (FileNotFoundException e) {
			return false;
		}
	}
}
