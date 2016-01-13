package application;

import java.util.ArrayList;

import javafx.scene.Node;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.RadioButton;
import javafx.scene.control.Slider;
import javafx.scene.control.Spinner;
import javafx.scene.control.TextArea;
import javafx.scene.control.Toggle;

// This class is a singleton
public class UpdateUserInterfaceSingleton {
	private static UpdateUserInterfaceSingleton instance = null;

	/*
	 * Singleton so that we only have one element updating and accessing the UI
	 * data.
	 */
	public static UpdateUserInterfaceSingleton getInstance() {
		if (instance == null) {
			instance = new UpdateUserInterfaceSingleton();
		}
		return instance;
	}

	/**
	 * Selects an element from the UI, if the element exists. Otherwise throws
	 * an error
	 * 
	 * @throws application.SceneNotYetReadyException
	 * 
	 * @throws NullPointerException
	 *
	 */
	public Node findById(String id) throws SceneNotYetReadyException {
		if (id == null)
			throw new NullPointerException();
		// CSS way of selecting
		return Main.lookupObjectById("#" + id);
	}

	/**
	 * Updates the values in the UI based on the values of the attributes being
	 * stored in Attributes.
	 */
	@SuppressWarnings("unchecked")
	public void updateDefaults() {
		System.out.print("Updating fields...");

		/*
		 * Iterate through our Attributes (The values stored by the
		 * attributesCollection hashmap)
		 */
		ArrayList<AttributeFactorySingleton.Attribute> array = AttributeFactorySingleton
				.getInstance().getAllAttributes();
		for (AttributeFactorySingleton.Attribute attribute : array) {
			// Get the correct node
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
			 * Type Slider
			 */
			if (currentData instanceof Slider) {
				if (!attribute.getValue().equals("NULL")) {
					((Slider) currentData).setValue(Double.valueOf(attribute
							.getValue()));
				}

				/*
				 * Type Spinner
				 */
			} else if (currentData instanceof Spinner) {
				if (!attribute.getValue().equals("NULL")) {
					Spinner<Integer> s = (Spinner<Integer>) currentData;
					s.getValueFactory().setValue(
							Integer.valueOf(attribute.getValue()));
				}

				/*
				 * Type Combo Box
				 */
			} else if (currentData instanceof ComboBox<?>) {
				((ComboBox<String>) currentData).setValue(attribute.getValue());

				/*
				 * Type Text Area
				 */
			} else if (currentData instanceof TextArea) {

				/*
				 * Some text fields naturally hold "NULL", but that is useless
				 * to the user. So we'll give the user "" if the actual value is
				 * NULL. It makes the UI look cleaner too because they don't
				 * have to look at a bunch of null fields.
				 */
				if (attribute.getValue().equals("NULL")) {
					((TextArea) currentData).setText("");
				} else {
					((TextArea) currentData).setText(attribute.getValue());
				}

				/*
				 * Type checkbox
				 */
			} else if (currentData instanceof CheckBox) {
				if (attribute.getValue().equals("")) {
					((CheckBox) currentData).setSelected(false);
				} else {
					((CheckBox) currentData).setSelected(true);
				}

				/*
				 * This is a special case, we use aAlgor to designate the
				 * alignment algorithm, which is transcribed to the UI via radio
				 * buttons. However, it is stored in the R file as a string.
				 */
			} else if (attribute.getName().equals("aAlgor")) {
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
