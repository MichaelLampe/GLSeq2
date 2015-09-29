package application;

import java.util.ArrayList;

import javafx.scene.Node;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.RadioButton;
import javafx.scene.control.Slider;
import javafx.scene.control.TextArea;
import javafx.scene.control.Toggle;

// This class is a singleton
public class UpdateUI {
	private static UpdateUI instance = null;
	private static Attributes attributes;

	protected UpdateUI() {
		attributes = Attributes.getInstance();
	}

	public static UpdateUI getInstance() {
		if (instance == null) {
			instance = new UpdateUI();
		}
		return instance;
	}

	public Node findById(String id) {
		if (id == null)
			throw new NullPointerException();
		// CSS way of selecting
		return Main.scene.lookup("#" + id);
	}

	@SuppressWarnings("unchecked")
	public void updateDefaults() {
		System.out.print("Updating fields...");
		// Create an array that can be iterates through from our HashMap
		ArrayList<Attribute> array = new ArrayList<Attribute>(attributes.attributesCollection.values());
		for (Attribute attribute : array) {
			// Get the correct node
			Node currentData = findById(attribute.getUiName());
			// Sliders we just want to set our saved value as a double.
			if (currentData instanceof Slider) {
				if (!attribute.getValue().equals("NULL")) {
					((Slider) currentData).setValue(Double.valueOf(attribute.getValue()));
				}
				// Combo boxes we just set the saved value directly.
			} else if (currentData instanceof ComboBox<?>) {
				((ComboBox<String>) currentData).setValue(attribute.getValue());
			} else if (currentData instanceof TextArea) {
				// Some text fields naturally hold "NULL", but that is useless
				// to the
				// user.
				// So we'll give the user "" if the actual value is NULL
				if (attribute.getValue().equals("NULL")) {
					((TextArea) currentData).setText("");
				} else {
					((TextArea) currentData).setText(attribute.getValue());
				}
				// Unchecked if empty string, otherwise check it
			} else if (currentData instanceof CheckBox) {
				if (attribute.getValue().equals("")) {
					((CheckBox) currentData).setSelected(false);
				} else {
					((CheckBox) currentData).setSelected(true);
				}
			} else if (attribute.getName().equals("aAlgor")) {
				// We go through all the radio buttons here to see which one
				// matches
				// The value of aAlgor and select that one.
				for (Toggle tog : MainPageController.group.getToggles()) {
					RadioButton button = (RadioButton) tog;
					if (button.getText().equals(attribute.getValue())) {
						button.setSelected(true);
					}
				}
			}
		}
		System.out.println("...Fields updated!");
	}
}
