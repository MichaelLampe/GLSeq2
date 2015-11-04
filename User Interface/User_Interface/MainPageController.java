package application;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.ResourceBundle;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import javafx.application.Platform;
import javafx.beans.binding.Bindings;
import javafx.beans.binding.BooleanBinding;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.ObservableList;
import javafx.concurrent.Task;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.event.EventType;
import javafx.fxml.Initializable;
import javafx.geometry.Insets;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonBar.ButtonData;
import javafx.scene.control.ButtonType;
import javafx.scene.control.CheckBox;
import javafx.scene.control.CheckBoxTreeItem;
import javafx.scene.control.Dialog;
import javafx.scene.control.Label;
import javafx.scene.control.PasswordField;
import javafx.scene.control.RadioButton;
import javafx.scene.control.SingleSelectionModel;
import javafx.scene.control.Tab;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.control.Toggle;
import javafx.scene.control.ToggleGroup;
import javafx.scene.control.Tooltip;
import javafx.scene.control.TreeCell;
import javafx.scene.control.TreeItem;
import javafx.scene.control.TreeView;
import javafx.scene.control.cell.CheckBoxTreeCell;
import javafx.scene.layout.GridPane;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.stage.FileChooser.ExtensionFilter;
import javafx.util.Callback;
import javafx.util.Pair;

public final class MainPageController extends MainPageItems implements
		Initializable {

	// Toggle groups to keep radio buttons together
	public static final ToggleGroup group = new ToggleGroup();
	private static final ToggleGroup referenceGroup = new ToggleGroup();
	private static final ToggleGroup fastqFiles = new ToggleGroup();
	private static final ToggleGroup samGroup = new ToggleGroup();

	// The head tree items of both of alignment and counting for advanced
	// options
	private static TreeItem<String> alignment;
	private static TreeItem<String> counting;
	private static TreeItem<String> data;

	// Hashes that contain all of the possible tree branches that could be added
	// to the map
	private static HashMap<String, TreeItem<String>> alignOptions = new HashMap<String, TreeItem<String>>();
	private static HashMap<String, TreeItem<String>> countOptions = new HashMap<String, TreeItem<String>>();

	// A single glow server request with username and password per UI open
	private static GlowRequest request = null;

	public final static ArrayList<String> liblistData = new ArrayList<String>();

	// CONSTANTS
	private final String LOGIN_GLOW = "Login to GLOW";
	private final String REQUEST_GLOW = "Request from GLOW";

	@Override
	public void initialize(URL arg0, ResourceBundle arg1) {
		addListeners();

		constructRadioButtons();

		setupComboBoxes();

		setupStartRun();

		addScrollbarListeners();

		progressBarListener();

		runBindings();
		menuListeners();

		appendTooltips();

		alignAndCountingOptions();

		fileSelectTree();

	}

	private void fillOptions(String option) {
		RadioButton[] tempAlignmentGroup = { BWA, Bowtie, Bowtie2, Cushaw,
				Cushaw_Gpu, TopHat, Rockhopper, star, hisat, BWA };
	}

	private void constructRadioButtons() {
		// Setup the alignment button group
		RadioButton[] tempAlignmentGroup = { BWA, Bowtie, Bowtie2, Cushaw,
				Cushaw_Gpu, TopHat, Rockhopper, star, hisat, BWA };
		List<RadioButton> alignmentGroup = Arrays.asList(tempAlignmentGroup);
		alignmentGroup.forEach((radioButton) -> {
			radioButton.setToggleGroup(group);
		});
		alignmentGroup.get(0).setSelected(true);
		toggleCountingDep(group);

		// Setup the type of run button group
		RadioButton[] tempTypeOfRun = { indic_man, from_glow, de_novo };
		List<RadioButton> typeOfRun = Arrays.asList(tempTypeOfRun);
		typeOfRun.forEach((radioButton) -> {
			radioButton.setToggleGroup(referenceGroup);
		});
		typeOfRun.get(0).setSelected(true);

		// Sam File Choice
		RadioButton[] tempSamSource = { man_sam, load_prev };
		List<RadioButton> samSource = Arrays.asList(tempSamSource);
		samSource.forEach((radioButton) -> {
			radioButton.setToggleGroup(samGroup);
		});
		samSource.get(0).setSelected(true);

		// Bind the GLOW button so that it is only enabled when GLOW is
		// selected.
		if (request == null) {
			ref_request_glow.setText(LOGIN_GLOW);
			data_request_glow.setText(LOGIN_GLOW);
		}

		ref_request_glow.disableProperty().bind(
				Bindings.when(from_glow.selectedProperty()).then(false)
						.otherwise(true));
		data_request_glow.disableProperty().bind(
				Bindings.when(load_glow.selectedProperty()).then(false)
						.otherwise(true));
		experiment_id.disableProperty().bind(
				Bindings.when(load_glow.selectedProperty()).then(false)
						.otherwise(true));

		// Setup where to get the data from
		RadioButton[] tempGlow = { man_fastq, load_glow };
		List<RadioButton> glow = Arrays.asList(tempGlow);
		glow.forEach((radioButton) -> {
			radioButton.setToggleGroup(fastqFiles);
		});
		glow.get(0).setSelected(true);

		// If you add more items here, you need to add more cases for
		// "requestFor"
		// in the lambda function below
		Button[] tempButtonAction = { ref_request_glow, data_request_glow };
		List<Button> buttonAction = Arrays.asList(tempButtonAction);
		buttonAction.parallelStream().forEach((button) -> {
			// Make sure this has a TextArea assigned to it!
				final TextArea requestFor;
				if (button.equals(ref_request_glow)) {
					requestFor = genom_ref_id;
				} else if (button.equals(data_request_glow)) {
					requestFor = experiment_id;
				} else {
					// Only would get here if more items are added above.
					requestFor = null;
				}

				button.setOnAction(new EventHandler<ActionEvent>() {
					@Override
					public void handle(ActionEvent arg0) {
						if (request == null) {
							// Log the user in and change the button name
							createLoginDialog();
						} else {
							Task<Object> refReq = new GlowReferenceRequest(
									request, requestFor.getText());
							// Run on a different thread so it doesn't lock up
							// the
							// UI.
							new Thread(refReq).start();
						}
					}
				});
			});
	}

	private void addListeners() {
		/*
		 * Creates a list of all the directories that need to be checked to see
		 * if they exist as the user types in the directory. Then iterates over
		 * the list in parallel and adds a listener.
		 */
		TextArea[] tempDirectoriesMayExist = { scriptDirectory, trimPath,
				fastqcPath, picardToolsPath, bwaPath, bam2wigPath,
				rockhopperPath, cushawGpuPath, topHatPath, cushawPath,
				cushawIndexPath, hisatPath, starPath };

		List<TextArea> directoriesMayExist = Arrays
				.asList(tempDirectoriesMayExist);

		directoriesMayExist.parallelStream().forEach((fileName) -> {
			fileName.textProperty().addListener(new ChangeListener<String>() {
				@Override
				public void changed(ObservableValue<? extends String> arg0,
						String arg1, String arg2) {
					File dest = new File(fileName.getText());
					// Red if it can't be written to.
					String color = "red";
					if (dest.canWrite()) {
						// If you can write to, makes it green.
						color = "green";
					}
					// Just a nice little border
					fileName.setStyle("-fx-border-color: " + color
							+ " ; -fx-border-width: 2px ;");
				}
			});
		});

		/*
		 * Creates a list of all directories that need to be checked if they are
		 * writable. Then iterates over the list in parallel and adds a
		 * listener.
		 */
		TextArea[] tempDirectoryMayBeWritable = { storageDestination,
				destinationDirectory };
		List<TextArea> directoryMayBeWritable = Arrays
				.asList(tempDirectoryMayBeWritable);

		directoryMayBeWritable.parallelStream().forEach((fileName) -> {
			fileName.textProperty().addListener(new ChangeListener<String>() {
				@Override
				public void changed(ObservableValue<? extends String> arg0,
						String arg1, String arg2) {
					File dest = new File(fileName.getText());
					// Red if it can't be written to.
					String color = "red";
					if (dest.canWrite()) {
						// If you can write to, makes it green.
						color = "green";
					}

					// Just a nice little border
					fileName.setStyle("-fx-border-color: " + color
							+ " ; -fx-border-width: 2px ;");
				}
			});
		});

	}

	private void fileSelectTree() {
		data = new CheckBoxTreeItem<String>("Retrieved Files");
		selectedDataFiles
				.setCellFactory(new Callback<TreeView<String>, TreeCell<String>>() {
					@Override
					public TreeCell<String> call(TreeView<String> param) {
						return new CheckBoxTreeCell<String>() {
							@Override
							public void updateItem(String item, boolean empty) {
								super.updateItem(item, empty);
								// If there is no information for the Cell, make
								// it
								// empty
								if (empty) {
									setGraphic(null);
									setText(null);
									// Otherwise if it's not representation as
									// an item
									// of the tree
									// is not a CheckBoxTreeItem, remove the
									// checkbox
									// item
								} else if (!(getTreeItem() instanceof CheckBoxTreeItem)) {
									setGraphic(null);
								}
							}
						};
					}
				});
		selectedDataFiles.setRoot(data);
		selectedDataFiles.setShowRoot(false);
	}

	private void alignAndCountingOptions() {
		constructOptionsMap();
		appendOptionsToTreeItems();
		/*
		 * Sets up the radio button and check mark listeners for when different
		 * options are selected.
		 */
		addAlignmentAndCountingListeners();
		/*
		 * The head of the alignment tree
		 */
		alignment = new TreeItem<String>("Alignment");
		alignmentOptions.setRoot(alignment);
		alignmentOptions
				.setCellFactory(new Callback<TreeView<String>, TreeCell<String>>() {
					@Override
					public TreeCell<String> call(TreeView<String> param) {
						return new CheckBoxTreeCell<String>() {
							@Override
							public void updateItem(String item, boolean empty) {
								super.updateItem(item, empty);
								if (getTreeItem() instanceof SpecialInputTreeItem){
									// TODO: Do special stuff here to make the correct inputs possible at this given cell.
								
								}
								// If there is no information for the Cell, make
								// it
								// empty
								if (empty) {
									setGraphic(null);
									setText(null);
									// Otherwise if it's not representation as
									// an item
									// of the tree
									// is not a CheckBoxTreeItem, remove the
									// checkbox
									// item
								} else if (!(getTreeItem() instanceof CheckBoxTreeItem)) {
									setGraphic(null);
								}
							}
						};
					}
				});
		// This is just the "alignment" option we see at the head of the field.
		alignmentOptions.setShowRoot(true);
		/*
		 * The head of the counting tree
		 */
		counting = new TreeItem<String>("Counting");
		countingOptions.setRoot(counting);
		countingOptions
				.setCellFactory(new Callback<TreeView<String>, TreeCell<String>>() {
					@Override
					public TreeCell<String> call(TreeView<String> param) {
						return new CheckBoxTreeCell<String>() {
							@Override
							public void updateItem(String item, boolean empty) {
								super.updateItem(item, empty);
								// If there is no information for the Cell, make
								// it
								// empty
								if (empty) {
									setGraphic(null);
									setText(null);
									// Otherwise if it's not representation as
									// an item
									// of the tree
									// is not a CheckBoxTreeItem, remove the
									// checkbox
									// item
								} else if (!(getTreeItem() instanceof CheckBoxTreeItem)) {
									setGraphic(null);
								}
							}
						};
					}
				});
		countingOptions.setShowRoot(true);
	}

	/*
	 * METHOD
	 * 
	 * Goal: Setup the alignment and counting selector listeners to properly
	 * append and remove alignment and counting options from the tree views.
	 * 
	 * Result: A dynamic and responsive tree view to edit advanced options.
	 */
	private void addAlignmentAndCountingListeners() {
		// Alignment RadioButton
		RadioButton[] tempModifyAlignOptions = { Bowtie, Bowtie2, Cushaw, BWA,
				Cushaw_Gpu, TopHat, star, hisat, Rockhopper };
		List<RadioButton> modifyAlignOptions = Arrays
				.asList(tempModifyAlignOptions);

		modifyAlignOptions.parallelStream().forEach((radioButton) -> {
			// Adds the options of that button into the HashMap
				fillOptions(radioButton.getText());
				// Then add listeners to everything
				radioButton.selectedProperty().addListener(
						new ChangeListener<Boolean>() {
							@Override
							public void changed(
									ObservableValue<? extends Boolean> arg0,
									Boolean arg1, Boolean arg2) {
								if (radioButton.isSelected()) {
									// Rockhopper is special because it is both
									// alignment
									// and counting.
									alignment.getChildren().clear();
									alignment.getChildren().add(
											alignOptions.get(radioButton
													.getText()));
									if (!radioButton.getText().equals(
											Bowtie.getText())) {
										if (!radioButton.getText().equals(
												Bowtie2.getText())) {
											alignment.getChildren().remove(
													alignOptions.get("RSEM"));
										}
									}
								} else {
									alignment.getChildren().remove(
											alignOptions.get(radioButton
													.getText()));
								}
							}
						});
			});

		// Counting CheckBoxes
		CheckBox[] tempModifyCountOptions = { RSEM, HTseq, FeatureCounts,
				Cufflinks };
		List<CheckBox> modifyCountOptions = Arrays
				.asList(tempModifyCountOptions);

		modifyCountOptions.parallelStream().forEach(
				(checkBox) -> {
					checkBox.selectedProperty().addListener(
							new ChangeListener<Boolean>() {
								@Override
								public void changed(
										ObservableValue<? extends Boolean> arg0,
										Boolean arg1, Boolean arg2) {
									if (checkBox.isSelected()) {
										counting.getChildren().add(
												countOptions.get(checkBox
														.getText()));
									} else {
										counting.getChildren().remove(
												countOptions.get(checkBox
														.getText()));
									}
								}
							});
				});
	}

	/*
	 * METHOD
	 * 
	 * Goal: Construct two Hash maps that hash the String keywords to TreeItems
	 * that contain all the important options of a given protocol
	 */
	private void constructOptionsMap() {
		// Aligners
		alignOptions.put("Bowtie", new TreeItem<String>("Bowtie"));
		alignOptions.put("Bowtie2", new TreeItem<String>("Bowtie2"));
		alignOptions.put("BWA", new TreeItem<String>("BWA"));
		alignOptions.put("Cushaw", new TreeItem<String>("Cushaw"));
		alignOptions.put("Cushaw_GPU", new TreeItem<String>("Cushaw_GPU"));
		alignOptions.put("HISAT", new TreeItem<String>("HISAT"));
		alignOptions.put("Rockhopper", new TreeItem<String>("Rockhopper"));
		alignOptions.put("STAR", new TreeItem<String>("STAR"));
		alignOptions.put("TopHat", new TreeItem<String>("TopHat"));
		// Counting
		countOptions.put("Cufflinks", new TreeItem<String>("Cufflinks"));
		countOptions
				.put("FeatureCounts", new TreeItem<String>("FeatureCounts"));
		countOptions.put("HTSeq", new TreeItem<String>("HTSeq"));
		countOptions.put("RSEM", new TreeItem<String>("RSEM"));
	}

	private void appendOptionsToTreeItems() {
		System.out.println("Started loading XML");
		TreeViewOptionsLoader load = new TreeViewOptionsLoader();
		for (String keys : alignOptions.keySet()) {
			try {
				alignOptions.get(keys).getChildren()
						.add(load.createTreeElementfromXml(keys));
			} catch (DuplicateElementsInXmlError e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	/*
	 * METHOD
	 * 
	 * Goal: Adds tool tips to a number of the fields.
	 */
	private void appendTooltips() {
		// Data tab tool tips

		// Common options
		Tooltip storageToolTip = new Tooltip(
				"This should be a parent directory to store the final result long-term (i.e. project folder on a bigdata volume); the actual directory will be a subfolder named after the run ID");
		storageDestination.setTooltip(storageToolTip);

		Tooltip destinationToolTip = new Tooltip(
				"This should be a directory on a fast-access volume (such as your home directory on the scarcity server system");
		destinationDirectory.setTooltip(destinationToolTip);

		// Start from Fastq files
		Tooltip dataToolTip = new Tooltip(
				"The directory containing FASTQ files (either compressed or uncompressed; please select the compression options in the \"Algorithms > Pre-processing\" accordingly)");
		dataDirectory.setTooltip(dataToolTip);

		// Continue from pre-aligned reads
		Tooltip countableToolTip = new Tooltip(
				"Start run from previously generated countable.sam files.");
		countableSamDir.setTooltip(countableToolTip);

		// Data characteristics
		Tooltip strainToolTip = new Tooltip("Strain name of the organism");
		strain.setTooltip(strainToolTip);

		Tooltip libTooltip = new Tooltip(
				"\"F\" if mRNA sequence corresponds to strandedness of the reads (Single End libraries) or the first read in the pair (Paired End libraries). \"R\" in the opposite case.");
		libStrands.setTooltip(libTooltip);

		Tooltip qScoreTooltip = new Tooltip(
				"Format of the score describing the confidence of base calling");
		qualityScores.setTooltip(qScoreTooltip);

		Tooltip zipToolTip = new Tooltip(
				"Original data is supplied in a compressed (.fastq.gz) format.  Uncompressed data will have .fastq extension (if not, uncompressing at the pre-processing stage is needed)");
		unzipped.setTooltip(zipToolTip);

		Tooltip pairedToolTip = new Tooltip(
				"If the reads were constructed in a single- or paired-ended manner.");
		pairedEnd.setTooltip(pairedToolTip);

		Tooltip presplitToolTip = new Tooltip(
				"For paired-end libraries, data is split into left- and right-read FASTQ files for each library (if not, splitting at the pre-processing stage may be needed).");
		presplit.setTooltip(presplitToolTip);

	}

	/*
	 * All scroll bars have listeners that allow for dynamic text update to
	 * their associated labels. Those are all created here.
	 */
	private void addScrollbarListeners() {
		// Number of Cores
		numberCores.valueProperty().addListener(new ChangeListener<Object>() {
			@Override
			public void changed(ObservableValue<?> arg0, Object arg1,
					Object arg2) {
				int cores = (int) numberCores.getValue();
				String message;
				if (cores > 1) {
					message = "You are now using " + cores + " cores";
				} else {
					message = "You are now using " + cores + " core";
				}
				cores_label.textProperty().setValue(message);
				// Update dataprep
				int streams = (int) numberStreamsDataPrep.getValue();
				if (streams > 1) {
					message = streams
							+ " computational streams during alignment. "
							+ streams * cores + " cores will be used in total.";
				} else {
					if (numberCores.getValue() > 1) {
						message = streams
								+ " computational stream during alignment. "
								+ streams * cores
								+ " cores will be used in total.";
					} else {
						message = streams
								+ " computational stream during alignment. "
								+ streams * cores
								+ " core will be used in total.";
					}
				}
				dataprep_streams_label.textProperty().setValue(message);
				// Update number
				streams = (int) numberStreams.getValue();
				if (streams > 1) {
					message = streams
							+ " computational streams during alignment. "
							+ streams * cores + " cores will be used in total.";
				} else {
					if (numberCores.getValue() > 1) {
						message = streams
								+ " computational stream during alignment. "
								+ streams * cores
								+ " cores will be used in total.";
					} else {
						message = streams
								+ " computational stream during alignment. "
								+ streams * cores
								+ " core will be used in total.";
					}
				}
				alignment_streams_label.textProperty().setValue(message);
			}
		});

		// Dataprep Streams
		numberStreamsDataPrep.valueProperty().addListener(
				new ChangeListener<Object>() {
					@Override
					public void changed(ObservableValue<?> arg0, Object arg1,
							Object arg2) {
						int streams = (int) numberStreamsDataPrep.getValue();
						String message;
						int cores = (int) numberCores.getValue();
						if (streams > 1) {
							message = streams
									+ " computational streams during data preparation. "
									+ streams * cores
									+ " cores will be used in total.";
						} else {
							if (numberCores.getValue() > 1) {
								message = streams
										+ " computational stream during data preparation. "
										+ streams * cores
										+ " cores will be used in total.";
							} else {
								message = streams
										+ " computational stream during data preparation. "
										+ streams * cores
										+ " core will be used in total.";
							}
						}
						dataprep_streams_label.textProperty().setValue(message);
					}
				});

		// Alignment Streams
		numberStreams.valueProperty().addListener(new ChangeListener<Object>() {
			@Override
			public void changed(ObservableValue<?> arg0, Object arg1,
					Object arg2) {
				int streams = (int) numberStreams.getValue();
				int cores = (int) numberCores.getValue();
				String message;
				if (streams > 1) {
					message = streams
							+ " computational streams during alignment. "
							+ streams * cores + " cores will be used in total.";
				} else {
					if (numberCores.getValue() > 1) {
						message = streams
								+ " computational stream during alignment. "
								+ streams * cores
								+ " cores will be used in total.";
					} else {
						message = streams
								+ " computational stream during alignment. "
								+ streams * cores
								+ " core will be used in total.";
					}
				}
				alignment_streams_label.textProperty().setValue(message);
			}
		});

		// Trim Head
		trimHead.valueProperty().addListener(new ChangeListener<Object>() {
			@Override
			public void changed(ObservableValue<?> arg0, Object arg1,
					Object arg2) {
				int trim = (int) trimHead.getValue();
				String message;
				if (trim > 1) {
					message = "Head trim, " + trim + " base";
				} else {
					message = "Head trim, " + trim + " bases";
				}
				trim_head_label.textProperty().setValue(message);
			}
		});

		// Min Read Length
		minTrim.valueProperty().addListener(new ChangeListener<Object>() {
			@Override
			public void changed(ObservableValue<?> arg0, Object arg1,
					Object arg2) {
				int read = (int) minTrim.getValue();
				String message;
				if (read > 1) {
					message = "Minimal trimmed read length of " + read
							+ " base";
				} else {
					message = "Minimal trimmed read length of " + read
							+ " bases";
				}
				min_read_length_label.textProperty().setValue(message);
			}
		});
	}

	/*
	 * Combo boxes need to be setup with presets. Those are setup here.
	 */
	private void setupComboBoxes() {
		// Trim Reads
		readTrim.getItems().addAll("TRUE", "FALSE");
		readTrim.setValue("TRUE");
		// Strandedness
		libStrands.getItems().addAll("R", "F", "NULL");
		libStrands.setValue("R");
		// Sequencing Platforms
		seqPlatform.getItems().addAll("illumina", "capillary", "ls454",
				"solid", "helicos", "iontorrent", "pacbio");
		seqPlatform.setValue("illumina");
		// Quality Score
		qualityScores.getItems().addAll("phred33", "phred64");
		qualityScores.setValue("phred33");
		// Zipped
		unzipped.getItems().addAll("TRUE", "FALSE");
		unzipped.setValue("FALSE");
		// Endedness
		pairedEnd.getItems().addAll("TRUE", "FALSE");
		pairedEnd.setValue("TRUE");
		// Presplit
		presplit.getItems().addAll("TRUE", "FALSE");
		presplit.setValue("TRUE");
	}

	/*
	 * The start button starts stuff. 1) Generate Attribute File 2) Launch the
	 * script.
	 */
	private void setupStartRun() {
		start_run.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent arg0) {

				if (load_glow.isSelected()) {
					String rList = "c(";
					for (int i = 0; i < liblistData.size(); i++) {
						if (i == liblistData.size() - 1) {
							rList += "\"" + liblistData.get(i) + "\")";
						} else {
							rList += "\"" + liblistData.get(i) + "\", ";
						}
					}
					Attributes.getInstance().attributesCollection.get(
							AttributesJSON.libList.field_name).setValue(rList);
				} else {
					Attributes.getInstance().attributesCollection.get(
							AttributesJSON.libList.field_name).setValue("NULL");
				}

				AttributeActionsUI action = new AttributeActionsUI();
				UpdateAttribute.getInstance().updateAttributes();
				String attributeFileLocation = null;
				try {
					if (run_name.getText().equals("")) {
						attributeFileLocation = action.writeAttributesFile(
								null, null);
					} else {
						attributeFileLocation = action.writeAttributesFile(
								run_name.getText(), null);
					}
					System.out.println(attributeFileLocation);
					action.saveConfigFile("AttributesConfig.txt");
				} catch (IOException e) {
					e.printStackTrace();
				}
				// Creates a run instance
				Run currentRun = constructRun();
				// Grabs the args from that instance
				// If true, run on HTcondor, else do not
				List<String> args = currentRun.constructArgs(
						htcondor_check.isSelected(), attributeFileLocation);
				if (args != null) {
					// Might add some logic in the future that starts an SSH on
					// windows/mac
					Task<Object> task = new RunWorker(args);
					// Run on a different thread so it doesn't lock up the UI.
					new Thread(task).start();
				} else {
					System.out
							.println("Please select something for the script to do.");
				}
			}
		});
	}

	public void runBindings() {
		start_run.disableProperty().bind(
				Bindings.when(
				// Requires a run name
						run_name.textProperty()
								// There is some text
								.isNotEqualTo("")
								.and(
								// 45 char names max
								// Just makes the table easier to
								// handle.
								textCountLimit(run_name, 45))
								.and(textAlphanumeric(run_name))
								.and(noExistRun(run_name))
								.and(noExistScript(run_name))
								.and(
								// One of the boxes must be check to do
								// something.
								alignment_check.selectedProperty().or(
										data_prep_check
												.selectedProperty()
												.or(counting_check
														.selectedProperty())
												.or(collect_check
														.selectedProperty()))))
						.then(false).otherwise(true));
	}

	private BooleanBinding textCountLimit(TextArea run_name, int char_count) {
		BooleanBinding binding = Bindings.createBooleanBinding(() -> (run_name
				.getText().length() <= char_count), run_name.textProperty());
		return binding;
	}

	private BooleanBinding textAlphanumeric(TextArea run_name) {
		BooleanBinding binding = Bindings.createBooleanBinding(() -> (!run_name
				.getText().matches("^.*[^a-zA-Z0-9].*$")), run_name
				.textProperty());
		return binding;
	}

	private BooleanBinding noExistRun(TextArea run_name) {
		// Update on run name change
		BooleanBinding bindingRun = Bindings.createBooleanBinding(
				() -> (!new File(scriptDirectory.getText() + run_name.getText()
						+ ".RData").exists()), run_name.textProperty());
		return bindingRun;
	}

	private BooleanBinding noExistScript(TextArea run_name) {
		// Update on script location change
		BooleanBinding bindingScript = Bindings.createBooleanBinding(
				() -> (!new File(scriptDirectory.getText() + run_name.getText()
						+ ".RData").exists()), scriptDirectory.textProperty());
		return bindingScript;
	}

	private Run constructRun() {
		// Update no longer in use.
		boolean data_prep = data_prep_check.isSelected();
		boolean alignment = alignment_check.isSelected();
		boolean counting = counting_check.isSelected();
		boolean collect = collect_check.isSelected();
		return new Run(data_prep, alignment, counting, collect,
				run_name.getText());
	}

	private void toggleCountingDep(ToggleGroup group) {
		RSEM.disableProperty().bind(
				Bindings.when(
						Bowtie.selectedProperty()
								.or(Bowtie2.selectedProperty())).then(false)
						.otherwise(true));
		group.selectedToggleProperty().addListener(
				new ChangeListener<Toggle>() {
					@Override
					public void changed(
							ObservableValue<? extends Toggle> observable,
							Toggle oldValue, Toggle newValue) {
						if (!Bowtie.isSelected() && !Bowtie2.isSelected()) {
							RSEM.setSelected(false);
						}
					}
				});
	}

	/*
	 * A bit of code to make the scrollbar on the bottom progress as the user
	 * clicks on the tabs in the menu. When they reach the run page it will be
	 * 100%.
	 */
	private void progressBarListener() {
		SingleSelectionModel<Tab> current_tab = tab_check.getSelectionModel();
		double percent_per_index = (100.0 / tab_check.getTabs().size());
		current_tab.selectedIndexProperty().addListener(
				new ChangeListener<Number>() {
					@Override
					public void changed(ObservableValue<? extends Number> arg0,
							Number arg1, Number arg2) {
						double currentPercent = ((current_tab
								.getSelectedIndex() + 2) * percent_per_index) / 100;
						footer_progress.setProgress(currentPercent);
					}
				});
	}

	private void logoutGlow() {
		request.logout();
		Main.stage
				.setTitle("|   GLSeq2 User Interface   -   Not logged into GLOW   |");
	}

	/*
	 * The menu on the upper part of the UI is setup here, so that it correctly
	 * interacts with the rest of the UI.
	 */
	private void menuListeners() {
		toggle_advanced.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				adv_opts.setDisable(!adv_opts.isDisabled());
			}
		});

		close_program.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				Platform.exit();
				System.exit(0);
			}
		});

		login_menu.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				if (request == null) {
					createLoginDialog();
				} else {
					logoutGlow();
					createLoginDialog();
				}
			}
		});

		logout_menu.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				// Only try to logout if currently logged in
				if (request != null) {
					logoutGlow();
				}
			}
		});
		// Only allows user to select R files
		ExtensionFilter R = new ExtensionFilter("R Attribute Files (*.R)",
				"*.R");
		open_att.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				FileChooser fc = new FileChooser();
				fc.getExtensionFilters().add(R);
				fc.setSelectedExtensionFilter(R);
				fc.setTitle("Open attribute file");
				File att_file = fc.showOpenDialog(Main.stage);
				if (att_file != null) {
					AttributeActions action = new AttributeActions();
					try {
						action.setAttributes(att_file);
					} catch (FileNotFoundException e) {
						System.out.println("Attribute File not found!");
					}
					String name = att_file.getName().split(".R")[0];
					run_name.setText(name);
					run_name.setDisable(true);
					// Remove binding and make sure the start button is active
					// Need to remove binding to directly set the button as not
					// disabled.
					start_run.disableProperty().unbind();
					start_run.setDisable(false);
				}
			}
		});

		template_att.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				FileChooser fc = new FileChooser();
				fc.getExtensionFilters().add(R);
				fc.setSelectedExtensionFilter(R);
				fc.setTitle("Open attribute file");
				File att_file = fc.showOpenDialog(Main.stage);
				try {
					AttributeActions action = new AttributeActions();
					action.setAttributes(att_file);
				} catch (FileNotFoundException e) {
					System.out.println("Attribute file not found!");
				}

			}
		});

		show_hide.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				ObservableList<Tab> tabs = tab_check.getTabs();
				if (tabs.contains(adv_opts)) {
					tabs.remove(adv_opts);
				} else {
					tabs.add(adv_opts);
				}
			}
		});

		save_curr.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				DirectoryChooser dc = new DirectoryChooser();
				dc.setTitle("Save attribute file");
				File directory = dc.showDialog(Main.stage);
				if (directory != null) {
					AttributeActions action = new AttributeActions();
					try {
						action.writeAttributesFile(null,
								directory.getAbsolutePath());
					} catch (IOException e) {
						// Should be fine with selector
					}
				}
			}
		});
	}

	/*
	 * Creates a login dialog
	 */
	private void createLoginDialog() {
		// Creates a new dialog which is just a popup
		Dialog<Pair<String, String>> login = new Dialog<>();
		ButtonType loginButtonType = new ButtonType("Login", ButtonData.OK_DONE);
		login.getDialogPane().getButtonTypes()
				.addAll(loginButtonType, ButtonType.CANCEL);
		login.setTitle("Login Dialog");
		login.setHeaderText("GLOW Login");

		// Layout from: http://code.makery.ch/blog/javafx-dialogs-official/
		GridPane grid = new GridPane();
		grid.setHgap(10);
		grid.setVgap(10);
		grid.setPadding(new Insets(20, 150, 10, 10));

		TextField username = new TextField();
		username.setPromptText("GLOW Username");
		PasswordField password = new PasswordField();
		password.setPromptText("GLOW Password");

		grid.add(new Label("GLOW Username:"), 0, 0);
		grid.add(username, 1, 0);
		grid.add(new Label("GLOW Password:"), 0, 1);
		grid.add(password, 1, 1);

		// Add everything to the dialog page
		login.getDialogPane().setContent(grid);

		// This returns the username and password as a pair so we can return
		// both
		// easily.
		login.setResultConverter(dialogButton -> {
			if (dialogButton == loginButtonType) {
				return new Pair<>(username.getText(), password.getText());
			}
			return null;
		});

		login.show();

		// The login button
		Button login_button = (Button) login.getDialogPane().lookupButton(
				loginButtonType);

		// Login listener
		login_button
				.addEventFilter(EventType.ROOT, e -> {
					if (e.getEventType().equals(ActionEvent.ACTION)) {
						e.consume();
						GlowRequest req = new GlowRequest(username.getText(),
								password.getText());
						if (req.requestCookie()) {

							// Login worked, assign it as the current login
							// session
						request = req;
						Task<Object> t = request;
						// Pings the server to verify connection
						new Thread(t).start();

						// Give a visual clue at the top of the UI that they are
						// currently
						// logged into the server. They are actually never
						// really
						// "logged in",
						// we just keep using their initial credentials to send
						// further
						// requests and automatically renew the credentials as
						// necessary.
						Main.stage
								.setTitle("|   GLSeq2 User Interface   -  Currently logged in as "
										+ username.getText() + "   |");

						/*
						 * If the user has logged in and not entered a
						 * destination directory, we'll just auto-populate this
						 * as being the base of their own directory
						 */
						data_request_glow.setText(REQUEST_GLOW);
						ref_request_glow.setText(REQUEST_GLOW);
						if (destinationDirectory.getText() == "") {
							destinationDirectory.setText("/home/GLBRCORG/"
									+ username.getText());
						}

						// Close the dialog as it's job is now complete
						login.close();

						// If the login failed
					} else {
						// http://stackoverflow.com/questions/29911552/how-to-shake-a-login-dialog-in-javafx-8
						try {
							/*
							 * Do a little shake animation. Also clears both the
							 * fields and adds a red border around them to
							 * further indicate that their credentials were
							 * wrong.
							 */
							ShakeTransition anim = new ShakeTransition(login
									.getDialogPane(), null);
							username.setText("");
							username.setStyle("-fx-border-color: red ; -fx-border-width: 1px ;");
							password.setText("");
							password.setStyle("-fx-border-color: red ; -fx-border-width: 1px ;");
							anim.playFromStart();
						} catch (Exception a) {
							// Do nothing
							// Exceptions could be caused because of the null
							// passed
							// into the
							// animation and this just removes the stack traces
						}
					}
				}
			}	);
	}
}