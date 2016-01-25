package application;

import javafx.scene.control.TreeItem;

// This extension just exists to be able to identify this later.
public class SpecialInputTreeItem<T> extends TreeItem<String> {
	private String input_type;
	private String default_value;
	private AdvancedOptionsTreeCell watcher;

	public SpecialInputTreeItem(String input_type, String default_value) {
		super();
		this.input_type = input_type;
		this.default_value = default_value;
	}

	public String getInputType() {
		return this.input_type;
	}

	public String getDefaultValue() {
		return this.default_value;
	}

	public String checkIfCorrectCell(Object testCellType) {
		return null;
	}

	public AdvancedOptionsTreeCell getWatchingTreeCell() {
		return watcher;
	}

	public void setWatchingTreeCell(AdvancedOptionsTreeCell watch) {
		this.watcher = watch;
	}
	
	public void setDefaultValue(String new_value){
		this.default_value = new_value;
	}
}
