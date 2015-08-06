package application;

import java.util.ArrayList;

import javafx.scene.Node;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.RadioButton;
import javafx.scene.control.Slider;
import javafx.scene.control.TextArea;

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
    ArrayList<Attribute> array = new ArrayList<Attribute>(attributes.attributesCollection.values());
    // Iterate through the attributes
    for (Attribute attribute : array) {
      // Grab the correct node in the fxml file by name
      Node currentData = findById(attribute.getUiName());
      // Sliders are doubles, so let's convert them to integers (No trailing
      // 0s) and then strings
      if (currentData instanceof Slider) {
        attribute
            .setValue(String.valueOf(Integer.valueOf((int) ((Slider) currentData).getValue())));
        // Combo boxes we just want the selected item.
      } else if (currentData instanceof ComboBox<?>) {
        attribute.setValue(String.valueOf(((ComboBox<?>) currentData).getSelectionModel()
            .getSelectedItem().toString()));
        // Just grab the text from a text area
      } else if (currentData instanceof TextArea) {
        attribute.setValue(((TextArea) currentData).getText());
        // Check boxes are only used currently for the counting and run
        // settings.
        // Run settings are not in the attribute file, so let's just format
        // this as we need it for the counting
      } else if (currentData instanceof CheckBox) {
        if (((CheckBox) currentData).isSelected()) {
          attribute.setValue(attribute.getName());
        } else {
          attribute.setValue("");
        }
        // Another special case, the radio buttons we need to grab the text from
        // instead of match it to a value.
      } else if (attribute.getName().equals("aAlgor")) {
        RadioButton selectedButton = (RadioButton) MainPageController.group.getSelectedToggle();
        attribute.setValue(selectedButton.getText());
      }
    }
  }
}
