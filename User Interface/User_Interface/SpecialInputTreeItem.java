package application;

import javafx.scene.control.TreeItem;

// This extension just exists to be able to identify this later.
public class SpecialInputTreeItem<T> extends TreeItem<String> {
	String input_type;
	public SpecialInputTreeItem(String input_type) {
		this.input_type = input_type;
	}

}
