package application;

import java.util.ArrayList;

import javafx.scene.Node;
import javafx.scene.control.CheckBox;
import javafx.scene.control.CheckBoxTreeItem;
import javafx.scene.control.ComboBox;
import javafx.scene.control.RadioButton;
import javafx.scene.control.Slider;
import javafx.scene.control.Spinner;
import javafx.scene.control.TextArea;
import javafx.scene.control.TreeItem;

// This class is a singleton, no need to have multiple attempts to update
public class UpdateAttributeSingleton {
	/*
	 * Singleton stuff.
	 */
	private static UpdateAttributeSingleton instance = null;

	public static UpdateAttributeSingleton getInstance() {
		if (instance == null) {
			instance = new UpdateAttributeSingleton();
		}
		return instance;
	}

	private Node findById(String id) throws SceneNotYetReadyException {
		if (id == null)
			throw new NullPointerException();

		/*
		 * The # looks for an ID, the id is a string value. Just regular CSS. A
		 * scene needs to exist for this to work.
		 */
		return Main.lookupObjectById("#" + id);
	}

	protected void updateAttributes() {
		// Iterate through the attributes
		for (AttributeFactorySingleton.Attribute attribute : AttributeFactorySingleton
				.getInstance().getAllAttributes()) {

			/*
			 * Grab the correct node in the Fxml file by name
			 */
			Node currentData;
			try {
				currentData = findById(attribute.getUiName());

				/*
				 * If the scene is not yet ready, no use trying to keep updating
				 * stuff.
				 */
			} catch (SceneNotYetReadyException e) {
				return;
			}

			/*
			 * Sliders are doubles, so let's convert them to integers (No
			 * trailing 0s) and then strings
			 */
			if (currentData instanceof Slider) {
				attribute.setValue(String.valueOf(Integer
						.valueOf((int) ((Slider) currentData).getValue())));

				/*
				 * Get spinner numeric values.
				 */
			} else if (currentData instanceof Spinner) {
				@SuppressWarnings("unchecked")
				Spinner<Integer> s = (Spinner<Integer>) currentData;
				attribute.setValue(String.valueOf(s.getValue()));

				/*
				 * Select the the selected item from teh combo box.
				 */
			} else if (currentData instanceof ComboBox<?>) {
				attribute.setValue(String.valueOf(((ComboBox<?>) currentData)
						.getSelectionModel().getSelectedItem().toString()));

				// Just grab the text from a text area
			} else if (currentData instanceof TextArea) {
				attribute.setValue(((TextArea) currentData).getText());

				/*
				 * Check boxes are only used currently for the counting and run
				 * settings. Run settings are not in the attribute file, so
				 * let's just format this as we need it for the counting
				 */
			} else if (currentData instanceof CheckBox) {
				if (((CheckBox) currentData).isSelected()) {
					attribute.setValue(attribute.getName());
				} else {
					attribute.setValue("");
				}

				/*
				 * Another special case, the radio buttons we need to grab the
				 * text from instead of match it to a value.
				 */
			} else if (attribute.getName().equals("aAlgor")) {
				RadioButton selectedButton = (RadioButton) MainPageController.group
						.getSelectedToggle();
				attribute.setValue(selectedButton.getText());
			}
		}
	}

	public void updateLiblist(ArrayList<String> liblistData) {
		if (liblistData.size() > 0) {
			String listOfFilesInR = "c(";
			for (int i = 0; i < liblistData.size(); i++) {
				if (i == liblistData.size() - 1) {
					listOfFilesInR += "\"" + liblistData.get(i) + "\")";
				} else {
					listOfFilesInR += "\"" + liblistData.get(i) + "\", ";
				}
			}
			try {
				AttributeFactorySingleton.getInstance().setAttributeValue(
						"libList", listOfFilesInR);
			} catch (NoSuchKeyInAttributeFactoryException e) {

				/*
				 * Print out error here as we can't throw from inside handler.
				 */
				e.printStackTrace();
			}
		} else {
			try {
				AttributeFactorySingleton.getInstance().setAttributeValue(
						"libList", "NULL");
			} catch (NoSuchKeyInAttributeFactoryException e) {

				/*
				 * Print out error here as we can't throw from inside handler.
				 */
				e.printStackTrace();
			}
		}
	}

	protected String constructSpecialArg(TreeItem<String> options,
			TreeViewOptionsLoader commandHolder) {

		/*
		 * String will be extended so use string builder.
		 */
		StringBuilder alignSpecialOptions = new StringBuilder("");

		/*
		 * If a check box is selected, add it to the strings so that it is
		 * considered an option.
		 */
		for (int i = 0; i < options.getChildren().size(); i++) {
			CheckBoxTreeItem<String> checkVariant = (CheckBoxTreeItem<String>) options
					.getChildren().get(i);
			if (checkVariant.isSelected()) {
				String optionName = options.getChildren().get(i).getValue();
				String command = TreeViewOptionsLoader.getCommandMap().get(
						optionName);
				alignSpecialOptions.append(command + " ");

				if (options.getChildren().get(i).getChildren().size() > 0) {

					/*
					 * Where the value is stored in the graphic. Either a
					 * spinner or a comboBox
					 */
					SpecialInputTreeItem<String> value = (SpecialInputTreeItem<String>) options
							.getChildren().get(i).getChildren().get(0);

					alignSpecialOptions.append(value.getDefaultValue() + " ");
				}
			}
		}

		/*
		 * Convert string builder to string upon return.
		 */
		return alignSpecialOptions.toString();
	}

	protected void addSpecialArgs(TreeItem<String> alignment,
			TreeItem<String> counting, TreeViewOptionsLoader commandHolder) {

		/*
		 * Alignment
		 */

		/*
		 * Root -> Aligner Name
		 */
		TreeItem<String> aligner = alignment.getChildren().get(0);

		/*
		 * Aligner Name -> Options
		 */
		TreeItem<String> options = aligner.getChildren().get(0);

		/*
		 * All the options
		 */
		try {
			AttributeFactorySingleton.getInstance().setAttributeValue(
					AttributeFactorySingleton.ALIGNMENT_SPECIAL_OPTIONS_NAME,
					constructSpecialArg(options, commandHolder));
		} catch (NoSuchKeyInAttributeFactoryException e1) {
			// Couldn't find the key.
		}

		/*
		 * Counting
		 */
		for (TreeItem<String> counter : counting.getChildren()) {
			TreeItem<String> option = counter.getChildren().get(0);
			/*
			 * Assign special options for counting protocols
			 */
			try {
				if (counter.getValue().equals("HTSeq")) {
					AttributeFactorySingleton
							.getInstance()
							.setAttributeValue(
									AttributeFactorySingleton.HTSEQ_SPECIAL_OPTIONS_NAME,
									constructSpecialArg(option, commandHolder));
				} else if (counter.getValue().equals("RSEM")) {
					AttributeFactorySingleton
							.getInstance()
							.setAttributeValue(
									AttributeFactorySingleton.RSEM_SPECIAL_OPTIONS_NAME,
									constructSpecialArg(option, commandHolder));
				} else if (counter.getValue().equals("FeatureCounts")) {
					AttributeFactorySingleton
							.getInstance()
							.setAttributeValue(
									AttributeFactorySingleton.FEATURE_COUNTS_SPECIAL_OPTIONS_NAME,
									constructSpecialArg(option, commandHolder));
				} else if (counter.getValue().equals("Cufflinks")) {
					AttributeFactorySingleton
							.getInstance()
							.setAttributeValue(
									AttributeFactorySingleton.CUFFLINKS_SPECIAL_OPTIONS_NAME,
									constructSpecialArg(option, commandHolder));
				}
			} catch (NoSuchKeyInAttributeFactoryException e) {
				continue;
			}
		}
	}

	public void updateAllAttributes(TreeItem<String> alignment,
			TreeItem<String> counting, TreeViewOptionsLoader commandHolder,
			ArrayList<String> libListData) {
		UpdateAttributeSingleton.getInstance().updateLiblist(libListData);
		UpdateAttributeSingleton.getInstance().addSpecialArgs(alignment,
				counting, commandHolder);
		UpdateAttributeSingleton.getInstance().updateAttributes();

	}
}
