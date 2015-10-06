package application;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import javafx.application.Platform;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.scene.control.CheckBoxTreeItem;
import javafx.scene.control.TreeItem;
import javafx.scene.control.TreeView;

public class GlowDataRequest extends GlowRequest {
	private String glow_id;

	public GlowDataRequest(GlowRequest connection, String glow_id) {
		super(connection.getUsername(), connection.getPassword());
		this.glow_id = glow_id;
		this.connected = connection.getConnected();
	}

	@SuppressWarnings("serial")
	private List<String> requestData(String glow_id) {
		ArrayList<String> command = new ArrayList<String>() {
			{
				add("curl");
				add("--cookie");
				add(cookie);
				add("--data");
				add("experiment_id=" + glow_id);
				add(glow_server + "get_experiment_data_by_id");
			}
		};
		return command;
	}

	@SuppressWarnings("unchecked")
	private ArrayList<CheckBoxTreeItem<String>> addFilesToUi(List<String> files) {
		try {
			TreeView<String> fileView = (TreeView<String>) UpdateUI.getInstance().findById("selectedDataFiles");
			for (String file : files) {
				int i = file.lastIndexOf("/");
				String filePath = file.substring(0, i + 1);
				String fileName = file.substring(i + 1);

				CheckBoxTreeItem<String> filePathBranch = null;
				CheckBoxTreeItem<String> fileNameBranch = new CheckBoxTreeItem<String>(fileName);
				for (TreeItem<String> path : fileView.getRoot().getChildren()) {
					if (path.getValue().equals(filePath)) {
						filePathBranch = (CheckBoxTreeItem<String>) path;
					}
				}

				// If still null, just makes a new one and adds it
				if (filePathBranch == null) {
					filePathBranch = new CheckBoxTreeItem<String>(filePath);
					fileView.getRoot().getChildren().add(filePathBranch);
				}

				// Add the file name to the path branch
				filePathBranch.getChildren().add(fileNameBranch);

				// Final reference for the listener
				final CheckBoxTreeItem<String> b = filePathBranch;

				// If we want to check the selected property we'll need to make
				// it into
				// a check box tree item

				b.selectedProperty().addListener(new ChangeListener<Boolean>() {
					@Override
					public void changed(ObservableValue<? extends Boolean> arg0, Boolean oldValue, Boolean newValue) {

						/*
						 * This gets activated when a change happens Changes all
						 * the TreeItems to be equal to their path branch
						 */

						for (TreeItem<String> val : b.getChildren()) {
							// Cast to correct type
							CheckBoxTreeItem<String> checks = (CheckBoxTreeItem<String>) val;
							// Set equal to the new value
							checks.setSelected(newValue);
						}
					}
				});
			}

			ArrayList<CheckBoxTreeItem<String>> fileBoxes = new ArrayList<CheckBoxTreeItem<String>>();
			for (TreeItem<String> path : fileView.getRoot().getChildren()) {
				for (TreeItem<String> file : path.getChildren()) {
					fileBoxes.add((CheckBoxTreeItem<String>) file);
				}
			}
			return fileBoxes;
		} catch (Exception e) {
			e.printStackTrace();
		}

		return null;
	}

	private void addBoxListeners(ArrayList<CheckBoxTreeItem<String>> files) {
		for (CheckBoxTreeItem<String> file : files) {
			file.selectedProperty().addListener(new ChangeListener<Boolean>() {
				@Override
				public void changed(ObservableValue<? extends Boolean> arg0, Boolean oldValue, Boolean newValue) {
					// Add and remove the string depending on if it is checked
					// or not
					if (newValue) {
						MainPageController.liblistData.add(file.getParent().getValue() + file.getValue());
					} else {
						MainPageController.liblistData.remove(file.getParent().getValue() + file.getValue());
					}
				}
			});
		}
	}

	@Override
	protected Object call() {
		if (renewCookieExpiration()) {
			Process process = null;
			try {
				process = new ProcessBuilder(requestData(glow_id)).start();
			} catch (IOException e2) {
				e2.printStackTrace();
			}
			String response = "";
			String line = null;
			try (BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
				while ((line = bufferedReader.readLine()) != null) {
					response += line;
				}
			} catch (IOException e2) {
				//
			}
			try {
				process.waitFor();
			} catch (InterruptedException e1) {
				//
			}
			ArrayList<String> foundFiles = new ArrayList<String>();
			String[] endingSplit = response.split("</file_path>");
			try {
				for (String phrase : endingSplit) {
					// Adds what was between the <file_path> tags
					String temp = phrase.split("<datafile")[1];
					try {
						String file = temp.split("<file_path>")[1];
						file = file.trim();
						if (file.endsWith(".gz")) {
							foundFiles.add(file);
						} else if (file.endsWith(".fastq")) {
							foundFiles.add(file);
						} else if (file.endsWith(".fq")) {
							foundFiles.add(file);
						}
					} catch (IndexOutOfBoundsException e) {

					}
				}
			} catch (IndexOutOfBoundsException e) {
				// Do nothing
			}
			ArrayList<CheckBoxTreeItem<String>> fileBoxes = addFilesToUi(foundFiles);
			addBoxListeners(fileBoxes);
			Platform.runLater(() -> {
				UpdateUI.getInstance().updateDefaults();
			});
		} else {
			System.out.println("Did not run because no connection to glow was established.");
		}
		return null;
	}
}