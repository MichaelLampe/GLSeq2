package application;

import java.util.ArrayList;

import javafx.scene.Node;
import javafx.scene.control.CheckBox;
import javafx.scene.control.CheckBoxTreeItem;
import javafx.scene.control.ComboBox;
import javafx.scene.control.RadioButton;
import javafx.scene.control.Slider;
import javafx.scene.control.TextArea;
import javafx.scene.control.TreeItem;

// This class is a singleton, no need to have multiple attempts to update
public class UpdateAttribute {
	private static UpdateAttribute instance = null;
	private static Attributes attributes;

	protected UpdateAttribute() {
		attributes = Attributes.getInstance();
	}

	public static UpdateAttribute getInstance() {
		if (instance == null) {
			instance = new UpdateAttribute();
		}
		return instance;
	}

	private Node findById(String id) {
		if (id == null)
			throw new NullPointerException();
		// The # looks for an ID, the id is a string value.
		// Just regular CSS
		return Main.scene.lookup("#" + id);
	}

	protected void updateAttributes() {
		// Create a new array from the HashMap so we can iterate through it
		ArrayList<Attribute> array = new ArrayList<Attribute>(
				attributes.attributesCollection.values());

		// Iterate through the attributes
		for (Attribute attribute : array) {
			// Grab the correct node in the Fxml file by name
			Node currentData = findById(attribute.getUiName());
			/*
			 * Sliders are doubles, so let's convert them to integers (No
			 * trailing 0s) and then strings
			 */
			if (currentData instanceof Slider) {
				attribute.setValue(String.valueOf(Integer
						.valueOf((int) ((Slider) currentData).getValue())));
				// Combo boxes we just want the selected item.
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

	protected String constructSpecialArg(TreeItem<String> options,
			TreeViewOptionsLoader commandHolder) {
		StringBuilder alignSpecialOptions = new StringBuilder("");
		for (int i = 0; i < options.getChildren().size(); i++) {
			CheckBoxTreeItem<String> checkVariant = (CheckBoxTreeItem<String>) options
					.getChildren().get(i);
			if (checkVariant.isSelected()) {
				String optionName = options.getChildren().get(i).getValue();
				String command = TreeViewOptionsLoader.getCommandMap().get(optionName);
				alignSpecialOptions.append(command);
				alignSpecialOptions.append(" ");

				if (options.getChildren().get(i).getChildren().size() > 0) {
					// Where the value is stored in the graphic.
					// Either a spinner or a comboBox

					SpecialInputTreeItem<String> value = (SpecialInputTreeItem<String>) options
							.getChildren().get(i).getChildren().get(0);

					alignSpecialOptions.append(value.getDefaultValue());
					alignSpecialOptions.append(" ");
				}
			}
		}
		return alignSpecialOptions.toString();
	}

	protected void addSpecialArgs(TreeItem<String> alignment,
			TreeItem<String> counting, TreeViewOptionsLoader commandHolder) {

		/*
		 * Alignment
		 */
		// Root -> Aligner Name
		TreeItem<String> aligner = alignment.getChildren().get(0);
		// Aligner Name -> Options
		TreeItem<String> options = aligner.getChildren().get(0);
		// All the options
		Attribute align = Attributes.getInstance().attributesCollection
				.get("alignmentSpecialOptions");
		align.setValue(constructSpecialArg(options, commandHolder));

		/*
		 * Counting
		 */
		// HtSeq Special Options
		for (TreeItem<String> counter : counting.getChildren()) {
			TreeItem<String> option = counter.getChildren().get(0);

			if (counter.getValue().equals("HTSeq")) {
				Attribute htseq = Attributes.getInstance().attributesCollection
						.get("HtSeqSpecialOptions");
				htseq.setValue(constructSpecialArg(option, commandHolder));

			}
			if (counter.getValue().equals("RSEM")) {
				Attribute rsem = Attributes.getInstance().attributesCollection
						.get("RsemSpecialOptions");
				rsem.setValue(constructSpecialArg(option, commandHolder));
			}
			if (counter.getValue().equals("FeatureCounts")) {
				Attribute featureCounts = Attributes.getInstance().attributesCollection
						.get("FeatureCountsSpecialOptions");
				featureCounts.setValue(constructSpecialArg(option,
						commandHolder));
			}
			if (counter.getValue().equals("Cufflinks")) {
				Attribute cufflinks = Attributes.getInstance().attributesCollection
						.get("CufflinksSpecialOptions");
				cufflinks.setValue(constructSpecialArg(option, commandHolder));
			}
		}
	}
}
