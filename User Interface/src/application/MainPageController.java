package application;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.ResourceBundle;

import com.jcabi.ssh.Shell;

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
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.geometry.Insets;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonBar.ButtonData;
import javafx.scene.control.ButtonType;
import javafx.scene.control.CheckBox;
import javafx.scene.control.CheckBoxTreeItem;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Dialog;
import javafx.scene.control.Label;
import javafx.scene.control.MenuItem;
import javafx.scene.control.PasswordField;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.RadioButton;
import javafx.scene.control.SingleSelectionModel;
import javafx.scene.control.Slider;
import javafx.scene.control.Spinner;
import javafx.scene.control.SpinnerValueFactory;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
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

public final class MainPageController implements Initializable {

	/*
	 * All of our FXML elements
	 */
	@FXML
	protected Slider gtfFeatureColumn;
	@FXML
	protected Label gtf_column_label;
	@FXML
	protected Slider fragMaxLength;
	@FXML
	protected Label max_fragment_length_label;
	@FXML
	protected Slider numberCores;
	@FXML
	protected Label cores_label;
	@FXML
	protected Spinner<Integer> trimHead;
	@FXML
	protected Spinner<Integer> minTrim;
	@FXML
	protected Spinner<Integer> three_prime_trim_slide;
	@FXML
	protected Spinner<Integer> min_avg_qual_slide;
	@FXML
	protected Slider ciMem;
	@FXML
	protected Label max_buffer_label;
	@FXML
	protected Slider numberStreamsDataPrep;
	@FXML
	protected Label dataprep_streams_label;
	@FXML
	protected Slider numberStreams;
	@FXML
	protected Label alignment_streams_label;
	@FXML
	protected ComboBox<String> readTrim;
	@FXML
	protected ComboBox<String> aAlgor;
	@FXML
	protected ComboBox<String> libStrands;
	@FXML
	protected ComboBox<String> seqPlatform;
	@FXML
	protected ComboBox<String> qualityScores;
	@FXML
	protected ComboBox<String> unzipped;
	@FXML
	protected ComboBox<String> pairedEnd;
	@FXML
	protected ComboBox<String> presplit;
	@FXML
	protected TextArea dataDirectory;
	@FXML
	protected TextArea destinationDirectory;
	@FXML
	protected TextArea countableSamDir;
	@FXML
	protected TextArea prevRunDirectory;
	@FXML
	protected TextArea prevRunName;
	@FXML
	protected TextArea referenceGenome;
	@FXML
	protected TextArea referencFasta;
	@FXML
	protected TextArea referenceGff;
	@FXML
	protected TextArea idAttr;
	@FXML
	protected TextArea artificialFasta;
	@FXML
	protected TextArea scriptDirectory;
	@FXML
	protected TextArea trimPath;
	@FXML
	protected TextArea fastqcPath;
	@FXML
	protected TextArea picardToolsPath;
	@FXML
	protected TextArea bwaPath;
	@FXML
	protected TextArea bam2wigPath;
	@FXML
	protected TextArea rockhopperPath;
	@FXML
	protected TextArea cushawGpuPath;
	@FXML
	protected TextArea topHatPath;
	@FXML
	protected TextArea cushawPath;
	@FXML
	protected TextArea cushawIndexPath;
	@FXML
	protected TextArea destDirTest;
	@FXML
	protected TextArea hisatPath;
	@FXML
	protected TextArea starPath;
	@FXML
	protected TextArea storageDestination;
	@FXML
	protected TextArea run_name;
	@FXML
	protected CheckBox data_prep_check;
	@FXML
	protected CheckBox alignment_check;
	@FXML
	protected CheckBox counting_check;
	@FXML
	protected CheckBox collect_check;
	@FXML
	protected CheckBox htcondor_check;
	@FXML
	protected CheckBox python_module_glseq_source;
	@FXML
	protected CheckBox RSEM;
	@FXML
	protected CheckBox HTseq;
	@FXML
	protected CheckBox FeatureCounts;
	@FXML
	protected CheckBox Cufflinks;
	@FXML
	protected Button start_run;
	@FXML
	protected RadioButton BWA;
	@FXML
	protected RadioButton Bowtie;
	@FXML
	protected RadioButton Bowtie2;
	@FXML
	protected RadioButton Cushaw;
	@FXML
	protected RadioButton Cushaw_Gpu;
	@FXML
	protected RadioButton TopHat;
	@FXML
	protected RadioButton Rockhopper;
	@FXML
	protected RadioButton star;
	@FXML
	protected RadioButton hisat;
	@FXML
	protected RadioButton load_prev;
	@FXML
	protected RadioButton man_sam;
	@FXML
	protected TabPane tab_check;
	@FXML
	protected ProgressBar footer_progress;
	@FXML
	protected MenuItem toggle_advanced;
	@FXML
	protected MenuItem close_program;
	@FXML
	protected MenuItem open_att;
	@FXML
	protected MenuItem template_att;
	@FXML
	protected MenuItem save_curr;
	@FXML
	protected MenuItem show_hide;
	@FXML
	protected Tab adv_opts;
	@FXML
	protected TextArea run_id;
	@FXML
	protected TreeView<String> alignmentOptions;
	@FXML
	protected TreeView<String> countingOptions;
	@FXML
	protected RadioButton from_glow;
	@FXML
	protected RadioButton indic_man;
	@FXML
	protected RadioButton de_novo;
	@FXML
	protected TextArea genom_ref_id;
	@FXML
	protected Button ref_request_glow;
	@FXML
	protected Button workingDirectoryPicker;
	@FXML
	protected Button destinationDirectoryPicker;
	@FXML
	protected Button data_request_glow;
	@FXML
	protected RadioButton load_glow;
	@FXML
	protected RadioButton man_fastq;
	@FXML
	protected TextArea experiment_id;
	@FXML
	protected TreeView<String> selectedDataFiles;
	@FXML
	protected MenuItem login_menu;
	@FXML
	protected MenuItem logout_menu;

	// Toggle groups to keep radio buttons together
	public static final ToggleGroup group = new ToggleGroup();
	private static final ToggleGroup referenceGroup = new ToggleGroup();
	private static final ToggleGroup fastqFiles = new ToggleGroup();
	private static final ToggleGroup samGroup = new ToggleGroup();

	boolean runningOnLinux = System.getProperty("os.name", "generic").toLowerCase(Locale.ENGLISH).indexOf("nux") >= 0;

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

	/*
	 * This is only public to handle the async callback in GlowDataRequest.
	 */
	public final static ArrayList<String> liblistData = new ArrayList<String>();

	private final static TreeViewOptionsLoader load = new TreeViewOptionsLoader();

	// CONSTANTS
	private final String LOGIN_GLOW = "Login to GLOW";
	private final String REQUEST_GLOW = "Request from GLOW";
	// Condor GLBRC "Scarcity-cm.glbrc.org" ip address
	protected final String condorAddress = "144.92.98.39";

	public static void logoutFromGlow() {
		if (request != null) {
			request.logout();
		}
	}

	public String getUsername() {
		// try to log the user in
		if (request == null) {
			createLoginDialog();
			try {
				return request.getUsername();
			} catch (NullPointerException e) {
				return null;
			}
		} else {
			return request.getUsername();
		}
	}

	public String getPassword() {
		if (request == null) {
			createLoginDialog();
			return null;
		}
		return request.getPassword();
	}

	@Override
	public void initialize(URL arg0, ResourceBundle arg1) {
		addTextAreaListeners();

		appendSliderListeners();

		progressBarListener();

		alignAndCountingOptions();

		constructRadioButtons();

		setupComboBoxes();

		setupStartRun();

		runBindings();

		menuListeners();

		appendTooltipsElements();

		fileSelectTree();

		instantiateUISpinnersWithValuesAndRanges();
	}

	/*
	 * This deals with settings up all the spinners. Spinners require
	 * IntegerSpinnerValueFactories to work, which also require an input range.
	 * 
	 * As of 1/25 there is also a bug fix I apply here so that if a user edits
	 * manually (Clicks and edits) the update will be tracked
	 */
	private void instantiateUISpinnersWithValuesAndRanges() {
		SpinnerValueFactory<Integer> trimHeadFactory = new SpinnerValueFactory.IntegerSpinnerValueFactory(0, 100);
		trimHead.setValueFactory(trimHeadFactory);
		trimHead.getValueFactory().setValue(12);

		/*
		 * Fixes a bug where JavaFX Spinner doesn't update when edited manually.
		 */
		trimHead.focusedProperty().addListener((observable, oldValue, newValue) -> {
			trimHead.getValueFactory().setValue(Integer.valueOf(trimHead.getEditor().getText()));
		});

		SpinnerValueFactory<Integer> minTrimFactory = new SpinnerValueFactory.IntegerSpinnerValueFactory(0, 100);
		minTrim.setValueFactory(minTrimFactory);
		minTrim.getValueFactory().setValue(36);
		minTrim.focusedProperty().addListener((observable, oldValue, newValue) -> {
			minTrim.getValueFactory().setValue(Integer.valueOf(minTrim.getEditor().getText()));
		});

		SpinnerValueFactory<Integer> threePrimeTrimFactory = new SpinnerValueFactory.IntegerSpinnerValueFactory(0, 100);
		three_prime_trim_slide.setValueFactory(threePrimeTrimFactory);
		three_prime_trim_slide.getValueFactory().setValue(3);
		three_prime_trim_slide.focusedProperty().addListener((observable, oldValue, newValue) -> {
			three_prime_trim_slide.getValueFactory()
					.setValue(Integer.valueOf(three_prime_trim_slide.getEditor().getText()));
		});

		SpinnerValueFactory<Integer> minAvgQualFactory = new SpinnerValueFactory.IntegerSpinnerValueFactory(0, 100);
		min_avg_qual_slide.setValueFactory(minAvgQualFactory);
		min_avg_qual_slide.getValueFactory().setValue(30);
		min_avg_qual_slide.focusedProperty().addListener((observable, oldValue, newValue) -> {
			min_avg_qual_slide.getValueFactory().setValue(Integer.valueOf(min_avg_qual_slide.getEditor().getText()));
		});
	}

	private void constructRadioButtons() {
		/*
		 * Alignment group button
		 */
		RadioButton[] tempAlignmentGroup = { BWA, Bowtie, Bowtie2, Cushaw, Cushaw_Gpu, TopHat, Rockhopper, star, hisat,
				BWA };
		List<RadioButton> alignmentGroup = Arrays.asList(tempAlignmentGroup);
		alignmentGroup.forEach((radioButton) -> {
			radioButton.setToggleGroup(group);
		});
		alignmentGroup.get(0).setSelected(true);
		toggleCountingDep(group);

		/*
		 * Setup the type of run button group
		 */
		RadioButton[] tempTypeOfRun = { indic_man, from_glow, de_novo };
		List<RadioButton> typeOfRun = Arrays.asList(tempTypeOfRun);
		typeOfRun.forEach((radioButton) -> {
			radioButton.setToggleGroup(referenceGroup);
		});
		typeOfRun.get(0).setSelected(true);

		/*
		 * Setup different SAM file choices
		 */
		RadioButton[] tempSamSource = { man_sam, load_prev };
		List<RadioButton> samSource = Arrays.asList(tempSamSource);
		samSource.forEach((radioButton) -> {
			radioButton.setToggleGroup(samGroup);
		});
		samSource.get(0).setSelected(true);

		/*
		 * A custom binding so that the request button responds to the radio
		 * button
		 */
		if (request == null) {
			ref_request_glow.setText(LOGIN_GLOW);
			data_request_glow.setText(LOGIN_GLOW);
		}
		ref_request_glow.disableProperty()
				.bind(Bindings.when(from_glow.selectedProperty()).then(false).otherwise(true));
		data_request_glow.disableProperty()
				.bind(Bindings.when(load_glow.selectedProperty()).then(false).otherwise(true));
		experiment_id.disableProperty().bind(Bindings.when(load_glow.selectedProperty()).then(false).otherwise(true));

		/*
		 * Setup where to get the data from
		 */
		RadioButton[] tempGlow = { man_fastq, load_glow };
		List<RadioButton> glow = Arrays.asList(tempGlow);
		glow.forEach((radioButton) -> {
			radioButton.setToggleGroup(fastqFiles);
		});
		glow.get(0).setSelected(true);

		workingDirectoryPicker.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent arg0) {
				DirectoryChooser fileChooser = new DirectoryChooser();
				fileChooser.setTitle("Working Directory Picker");
				File selectedDirectory = fileChooser.showDialog(null);
				if (selectedDirectory != null) {
					destinationDirectory.setText(selectedDirectory.getAbsolutePath());
				}
			}
		});

		destinationDirectoryPicker.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent arg0) {
				DirectoryChooser fileChooser = new DirectoryChooser();
				fileChooser.setTitle("Destination Directory Picker");
				File selectedDirectory = fileChooser.showDialog(null);
				if (selectedDirectory != null) {
					storageDestination.setText(selectedDirectory.getAbsolutePath());
				}
			}
		});

		/*
		 * If you add more items here, you need to add more cases for
		 * "requestFor" in the lambda function below
		 */
		Button[] tempButtonAction = { ref_request_glow, data_request_glow };
		List<Button> buttonAction = Arrays.asList(tempButtonAction);

		// Add more cases here
		buttonAction.parallelStream().forEach((button) -> {
			// Make sure this has a TextArea assigned to it!
			final TextArea requestFor;
			if (button.equals(ref_request_glow)) {
				requestFor = genom_ref_id;
				button.setOnAction(new EventHandler<ActionEvent>() {
					@Override
					public void handle(ActionEvent arg0) {
						if (request == null) {
							// Log the user in and change the button name
							createLoginDialog();
						}
						// Starts right away if they actually did login.
						if (request != null) {
							Task<Object> refReq = new GlowReferenceRequest(request, requestFor.getText());
							/*
							 * Run on a different thread so it doesn't lock up
							 * the UI.
							 */
							new Thread(refReq).start();
						}
					}
				});
			} else if (button.equals(data_request_glow)) {
				requestFor = experiment_id;
				button.setOnAction(new EventHandler<ActionEvent>() {
					@Override
					public void handle(ActionEvent arg0) {
						if (request == null) {
							// Log the user in and change the button name
							createLoginDialog();
						}
						// Starts right away if they actually did login.
						if (request != null) {
							Task<Object> dataReq = new GlowDataRequest(request, requestFor.getText());
							/*
							 * Run on a different thread so it doesn't lock up
							 * the UI.
							 */
							new Thread(dataReq).start();
						}
					}
				});
			} else {
				// Only would get here if more items are added above.
				requestFor = null;
			}
		});
	}

	private void addTextAreaListeners() {
		/*
		 * Creates a list of all the directories that need to be checked to see
		 * if they exist as the user types in the directory. Then iterates over
		 * the list in parallel and adds a listener.
		 */
		TextArea[] tempDirectoriesMayExist = { scriptDirectory, trimPath, fastqcPath, picardToolsPath, bwaPath,
				bam2wigPath, rockhopperPath, cushawGpuPath, topHatPath, cushawPath, cushawIndexPath, hisatPath,
				starPath };

		List<TextArea> directoriesMayExist = Arrays.asList(tempDirectoriesMayExist);

		directoriesMayExist.parallelStream().forEach((fileName) -> {
			fileName.textProperty().addListener(new ChangeListener<String>() {
				@Override
				public void changed(ObservableValue<? extends String> arg0, String arg1, String new_file) {
					if (runningOnLinux) {
						File dest = new File(new_file);
						String color = "red";
						if (dest.isDirectory() || dest.exists()) {
							color = "green";
						}
						// Just a nice little border
						fileName.setStyle("-fx-border-color: " + color + " ; -fx-border-width: 2px ;");
					}
				}
			});
		});

		/*
		 * Adds the files in the given directory to the select file.
		 */
		dataDirectory.textProperty().addListener(new ChangeListener<String>() {

			@Override
			public void changed(ObservableValue<? extends String> arg0, String oldValue, String newValue) {
				File folder = new File(newValue);
				File[] files = folder.listFiles();

				// Ensure file is not null
				if (files == null) {
					return;
				}

				// Add directory to the head.
				String directory = folder.getAbsolutePath();
				TreeItem<String> itemDirectory = new TreeItem<String>(directory);

				// Add each file to the directory.
				for (File file : files) {
					if (file.isFile()) {
						if (file.getName().endsWith(".fq") || file.getName().endsWith(".fastq")
								|| file.getName().endsWith(".gz")) {
							CheckBoxTreeItem<String> checkbox = new CheckBoxTreeItem<String>(file.getName());
							checkbox.selectedProperty().addListener(new ChangeListener<Boolean>() {
								@Override
								public void changed(ObservableValue<? extends Boolean> arg0, Boolean oldValue,
										Boolean newValue) {
									/*
									 * Add and remove the string depending on if
									 * it is checked or not
									 */
									String directoryValue = checkbox.getParent().getValue();
									if (!directoryValue.endsWith("/")) {
										directoryValue += "/";
									}
									if (newValue) {
										MainPageController.liblistData.add(directoryValue + checkbox.getValue());
									} else {
										MainPageController.liblistData.remove(directoryValue

												+ checkbox.getValue());
									}
								}
							});
							itemDirectory.getChildren().add(checkbox);
						}
					}
				}
				// Reset the file select tree as we are going to add a new set
				// of items.
				fileSelectTree();
				selectedDataFiles.getRoot().getChildren().add(itemDirectory);
			}
		});

		/*
		 * Creates a list of all directories that need to be checked if they are
		 * writable. Then iterates over the list in parallel and adds a
		 * listener.
		 */
		TextArea[] tempDirectoryMayBeWritable = { storageDestination, destinationDirectory };
		List<TextArea> directoryMayBeWritable = Arrays.asList(tempDirectoryMayBeWritable);

		directoryMayBeWritable.parallelStream().forEach((fileName) ->

		{
			fileName.textProperty().addListener(new ChangeListener<String>() {
				@Override
				public void changed(ObservableValue<? extends String> arg0, String arg1, String arg2) {
					if (runningOnLinux) {
						File dest = new File(fileName.getText());
						String color = "red";
						if (dest.canWrite()) {
							color = "green";
						}

						// Just a nice little border
						fileName.setStyle("-fx-border-color: " + color + " ; -fx-border-width: 2px ;");
					}
				}
			});
		});

	}

	private void fileSelectTree() {
		// Remove all the items if there are any inside.
		for (String value : liblistData) {
			liblistData.remove(value);
		}
		data = new CheckBoxTreeItem<String>("Retrieved Files");
		selectedDataFiles.setCellFactory(new Callback<TreeView<String>, TreeCell<String>>() {
			@Override
			public TreeCell<String> call(TreeView<String> param) {
				return new CheckBoxTreeCell<String>() {
					@Override
					public void updateItem(String item, boolean empty) {
						super.updateItem(item, empty);
						/*
						 * If there is no information for the Cell, make it
						 * empty
						 */
						if (empty) {
							setGraphic(null);
							setText(null);
							/*
							 * Otherwise if it's not representation as an item
							 * of the tree is not a CheckBoxTreeItem, remove the
							 * checkbox item
							 */
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
		/*
		 * Sets up the radio button and check mark listeners for when different
		 * options are selected.
		 */
		addAlignmentAndCountingListeners();
		constructAlignmentAndCountingOptionsTreeItemMap();
		appendOptionsToTreeItems();
		/*
		 * The head of the alignment tree
		 */
		alignment = new TreeItem<String>("Alignment");
		alignmentOptions.setRoot(alignment);

		alignmentOptions.setCellFactory(new Callback<TreeView<String>, TreeCell<String>>() {
			@Override
			public TreeCell<String> call(TreeView<String> param) {
				return new AdvancedOptionsTreeCell();
			}
		});
		// This is just the "alignment" option we see at the head of the field.
		alignmentOptions.setShowRoot(true);
		/*
		 * The head of the counting tree
		 */
		counting = new TreeItem<String>("Counting");
		countingOptions.setRoot(counting);
		countingOptions.setCellFactory(new Callback<TreeView<String>, TreeCell<String>>() {
			@Override
			public TreeCell<String> call(TreeView<String> param) {
				return new AdvancedOptionsTreeCell();
			}
		});
		countingOptions.setShowRoot(true);
	}

	/*
	 * Adds various listeners to the alignment and coutning options so they are
	 * responsive to what is actually available.
	 */
	private void addAlignmentAndCountingListeners() {
		// Alignment RadioButton
		RadioButton[] tempModifyAlignOptions = { Bowtie, Bowtie2, Cushaw, BWA, Cushaw_Gpu, TopHat, star, hisat,
				Rockhopper };
		List<RadioButton> modifyAlignOptions = Arrays.asList(tempModifyAlignOptions);

		modifyAlignOptions.parallelStream().forEach((radioButton) -> {
			// Then add listeners to everything
			radioButton.selectedProperty().addListener(new ChangeListener<Boolean>() {
				@Override
				public void changed(ObservableValue<? extends Boolean> arg0, Boolean arg1, Boolean arg2) {
					if (radioButton.isSelected()) {
						// Rockhopper is special because it is both
						// alignment
						// and counting.
						alignment.getChildren().clear();
						alignment.getChildren().add(alignOptions.get(radioButton.getText()));
						if (!radioButton.getText().equals(Bowtie.getText())) {
							if (!radioButton.getText().equals(Bowtie2.getText())) {
								alignment.getChildren().remove(alignOptions.get("RSEM"));
							}
						}
					} else {
						alignment.getChildren().remove(alignOptions.get(radioButton.getText()));
					}
				}
			});
		});

		// Counting CheckBoxes
		CheckBox[] tempModifyCountOptions = { RSEM, HTseq, FeatureCounts, Cufflinks };
		List<CheckBox> modifyCountOptions = Arrays.asList(tempModifyCountOptions);
		modifyCountOptions.parallelStream().forEach((checkBox) -> {
			checkBox.selectedProperty().addListener(new ChangeListener<Boolean>() {

				@Override
				public void changed(ObservableValue<? extends Boolean> arg0, Boolean arg1, Boolean arg2) {
					String countOptionName = checkBox.getText();
					// Feature counts has a space, so remove the
					// space
					countOptionName = countOptionName.replace(" ", "");
					if (checkBox.isSelected()) {
						counting.getChildren().add(countOptions.get(countOptionName));
					} else {
						counting.getChildren().remove(countOptions.get(countOptionName));
					}
				}

			});
		});
	}

	/*
	 * Codes in the hash of all the options
	 */
	private void constructAlignmentAndCountingOptionsTreeItemMap() {
		/*
		 * Aligners
		 */
		alignOptions.put("Bowtie", new TreeItem<String>("Bowtie"));
		alignOptions.put("Bowtie2", new TreeItem<String>("Bowtie2"));
		alignOptions.put("BWA", new TreeItem<String>("BWA"));
		alignOptions.put("Cushaw", new TreeItem<String>("Cushaw"));
		alignOptions.put("Cushaw_GPU", new TreeItem<String>("Cushaw_GPU"));
		alignOptions.put("HISAT", new TreeItem<String>("HISAT"));
		alignOptions.put("Rockhopper", new TreeItem<String>("Rockhopper"));
		alignOptions.put("STAR", new TreeItem<String>("STAR"));
		alignOptions.put("TopHat", new TreeItem<String>("TopHat"));

		/*
		 * Counting
		 */
		countOptions.put("Cufflinks", new TreeItem<String>("Cufflinks"));
		countOptions.put("FeatureCounts", new TreeItem<String>("FeatureCounts"));
		countOptions.put("HTSeq", new TreeItem<String>("HTSeq"));
		countOptions.put("RSEM", new TreeItem<String>("RSEM"));
	}

	private void appendOptionsToTreeItems() {
		for (String keys : alignOptions.keySet()) {
			try {
				alignOptions.get(keys).getChildren().add(load.createTreeElementfromXml(keys));
			} catch (DuplicateElementsInXmlException e) {
				e.printStackTrace();
			}
		}

		for (String keys : countOptions.keySet()) {
			try {
				countOptions.get(keys).getChildren().add(load.createTreeElementfromXml(keys));
			} catch (DuplicateElementsInXmlException e) {
				e.printStackTrace();
			}
		}
	}

	/*
	
	 */
	private void appendTooltipsElements() {
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
		Tooltip countableToolTip = new Tooltip("Start run from previously generated countable.sam files.");
		countableSamDir.setTooltip(countableToolTip);

		Tooltip libTooltip = new Tooltip(
				"\"F\" if mRNA sequence corresponds to strandedness of the reads (Single End libraries) or the first read in the pair (Paired End libraries). \"R\" in the opposite case.");
		libStrands.setTooltip(libTooltip);

		Tooltip qScoreTooltip = new Tooltip("Format of the score describing the confidence of base calling");
		qualityScores.setTooltip(qScoreTooltip);

		Tooltip zipToolTip = new Tooltip(
				"Original data is supplied in a compressed (.fastq.gz) format.  Uncompressed data will have .fastq extension (if not, uncompressing at the pre-processing stage is needed)");
		unzipped.setTooltip(zipToolTip);

		Tooltip pairedToolTip = new Tooltip("If the reads were constructed in a single- or paired-ended manner.");
		pairedEnd.setTooltip(pairedToolTip);

		Tooltip presplitToolTip = new Tooltip(
				"For paired-end libraries, data is split into left- and right-read FASTQ files for each library (if not, splitting at the pre-processing stage may be needed).");
		presplit.setTooltip(presplitToolTip);
	}

	private void updateSliderFields() {
		int cores = (int) numberCores.getValue();
		int dataPrepStreams = (int) numberStreamsDataPrep.getValue();
		int alignmentStreams = (int) numberStreams.getValue();

		/*
		 * Cores message
		 */
		String coreMessage = "";
		if (cores > 1) {
			coreMessage = cores + " cores used for a single sample processing";
		} else {
			coreMessage = cores + " core used for a single sample processing";
		}
		cores_label.textProperty().setValue(coreMessage);

		/*
		 * Dataprep stream message
		 */
		String dataprepStreamMessage = "";
		if (dataPrepStreams > 1) {
			dataprepStreamMessage = dataPrepStreams + " samples processed in parallel during data preparation. "
					+ dataPrepStreams * cores + " cores will be used in total.";
		} else {
			if (numberCores.getValue() > 1) {
				dataprepStreamMessage = dataPrepStreams + " sample processed in parallel during data preparation. "
						+ dataPrepStreams * cores + " cores will be used in total.";
			} else {
				dataprepStreamMessage = dataPrepStreams + " sample processed in parallel during data preparation. "
						+ dataPrepStreams * cores + " core will be used in total.";
			}
		}
		dataprep_streams_label.textProperty().setValue(dataprepStreamMessage);

		String alignmentStreamMessage = "";
		if (alignmentStreams > 1) {
			alignmentStreamMessage = alignmentStreams + " samples processed in parallel during alignment. "
					+ alignmentStreams * cores + " cores will be used in total.";
		} else {
			if (numberCores.getValue() > 1) {
				alignmentStreamMessage = alignmentStreams + " sample processed in parallel during alignment. "
						+ alignmentStreams * cores + " cores will be used in total.";
			} else {
				alignmentStreamMessage = alignmentStreams + " sample processed in parallel during alignment. "
						+ alignmentStreams * cores + " core will be used in total.";
			}
		}
		alignment_streams_label.textProperty().setValue(alignmentStreamMessage);
	}

	/*
	 * All scroll bars have listeners that allow for dynamic text update to
	 * their associated labels. Those are all created here.
	 */
	private void appendSliderListeners() {
		// Number of Cores
		numberCores.valueProperty().addListener(new ChangeListener<Object>() {
			@Override
			public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
				updateSliderFields();
			}
		});

		// Dataprep Streams
		numberStreamsDataPrep.valueProperty().addListener(new ChangeListener<Object>() {
			@Override
			public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
				updateSliderFields();
			}
		});

		// Alignment Streams
		numberStreams.valueProperty().addListener(new ChangeListener<Object>() {
			@Override
			public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
				updateSliderFields();
			}
		});
	}

	/*
	 * Combo boxes need to be setup with presets. Those are setup here.
	 */
	private void setupComboBoxes() {
		// Trim Reads
		readTrim.getItems().addAll("TRUE", "FALSE");
		readTrim.setValue("FALSE");
		artificialFasta.disableProperty().set(true);
		minTrim.disableProperty().set(true);
		three_prime_trim_slide.disableProperty().set(true);
		trimHead.disableProperty().set(true);
		min_avg_qual_slide.disableProperty().set(true);
		readTrim.valueProperty().addListener(new ChangeListener<String>() {
			@Override
			public void changed(ObservableValue<? extends String> ov, String oldValue, String newValue) {
				if (newValue.equals("TRUE")) {
					artificialFasta.disableProperty().set(false);
					minTrim.disableProperty().set(false);
					three_prime_trim_slide.disableProperty().set(false);
					trimHead.disableProperty().set(false);
					min_avg_qual_slide.disableProperty().set(false);
				} else {
					artificialFasta.disableProperty().set(true);
					minTrim.disableProperty().set(true);
					three_prime_trim_slide.disableProperty().set(true);
					trimHead.disableProperty().set(true);
					min_avg_qual_slide.disableProperty().set(true);
				}
			}
		});
		// Strandedness
		libStrands.getItems().addAll("R", "F", "NULL");
		libStrands.setValue("R");
		// Sequencing Platforms
		seqPlatform.getItems().addAll("illumina", "capillary", "ls454", "solid", "helicos", "iontorrent", "pacbio");
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
				if (!runningOnLinux) {
					if (request == null) {
						createLoginDialog();
						// They didn't login so don't try to run
						if (request == null) {
							return;
						}
					}
				}
				/*
				 * Update absolutely everything.
				 */
				UpdateAttributeSingleton.getInstance().updateAllAttributes(alignment, counting, load, liblistData);
				if (python_module_glseq_source.isSelected()){
					try {
						UpdateAttributeSingleton.getInstance().clearBaseDir();
					} catch (NoSuchKeyInAttributeFactoryException e) {
						System.out.println("Base.dir not found");
					}
				}
				
				AttributeFileWriter attributeFileWriter = new AttributeFileWriter();

				/*
				 * Writes the attribute file in
				 */
				String attributeFileLocation = null;
				try {
					if (run_name.getText().equals("")) {
						attributeFileLocation = attributeFileWriter.writeAttributesFile(null, null);
					} else {
						attributeFileLocation = attributeFileWriter.writeAttributesFile(run_name.getText(), null);
					}
				} catch (IOException e) {
					e.printStackTrace();
				}

				/*
				 * Saves config file
				 */
				AttributeAndConfigFileHandler configFileSaver = new AttributeAndConfigFileHandler();
				try {
					configFileSaver.saveConfigFile("AttributesConfig.txt");
				} catch (IOException e) {
					e.printStackTrace();
				}

				/*
				 * Running on Linux
				 */
				if (runningOnLinux) {
					/*
					 * Creates a run instance
					 */
					RunInstantiator currentRun = constructRun();
					currentRun.constructArgs(htcondor_check.isSelected(), python_module_glseq_source.isSelected(), attributeFileLocation,
							scriptDirectory.getText());
					currentRun.start();

					String command = currentRun.returnArgString();
					Alert alert = new Alert(AlertType.INFORMATION);
					alert.setTitle("Job started");
					alert.setHeaderText("You have launched a new local linux pipeline job");
					alert.setContentText("Job launched with command \n'" + command + "'");
					alert.showAndWait();
				} else {
					/*
					 * Creates a run instance, condor always activated
					 */
					String dest = destinationDirectory.getText();

					if (!dest.endsWith("/")) {
						dest += "/";
					}

					String attributeFileLocationOnLinux = dest + run_name.getText() + ".R";
					;

					RunInstantiator currentRun = constructRun();
					currentRun.constructArgs(true, python_module_glseq_source.isSelected(), attributeFileLocationOnLinux, scriptDirectory.getText());
					String command = currentRun.returnArgString();
					// This throws a null pointer when you first hit it. Still
					// working on how I want to resolve it. The application
					// thread catches the pointer and doesn't crash and then the
					// user just needs to hit it again after login and it is all
					// ok, so it isn't that pressing of an issue.
					try {
						// Transfer attribute file to server
						System.out.println("Splitting command");
						String[] attributeFileLines = attributeFileWriter.fullFileAsString().split("\n");
						System.out.println("Finished splitting");

						// Send as one large formatted blob.
						StringBuilder finalSend = new StringBuilder();
						for (String line : attributeFileLines) {
							// Don't send comments
							if (!line.contains("#")) {
								// Append to file
								line += "\n";
								finalSend.append(line);
							}
						}
						Alert alertInitial = new Alert(AlertType.INFORMATION);
						alertInitial.setTitle("Sending job to server");
						alertInitial.setHeaderText("Transfering attribute file to server and instantiating job");
						alertInitial.setContentText(
								"Another alert will occur when job has been successfully submitted to server");
						// Nonblocking
						alertInitial.show();
						// Make sure the directory for the attribute file is there.
						new Thread(new Shell.Plain(request.establishSsh()).exec("mkdir -p " + dest)).start();
						
						String toFile = "echo '" + finalSend.toString() + "' >> " + attributeFileLocationOnLinux;
						
						System.out.println("Sent to the server: ");
						System.out.println(toFile);
						String response = new Shell.Plain(request.establishSsh()).exec(toFile);
						
						System.out.println("Server responded with: ");
						System.out.println(response);
						System.out.println("Finished writing attribute file");					
						String re = new Shell.Plain(request.establishSsh()).exec(command);//).start();
						System.out.println(re);
						System.out.println("Started job");
						Alert alertFinal = new Alert(AlertType.INFORMATION);
						alertFinal.setTitle("Job started");
						alertFinal.setHeaderText("You have launched a Condor job");
						alertFinal.setContentText("Job launched with command \n'" + command + "'");
						alertFinal.showAndWait();

					} catch (IOException e1) {
						// TODO Auto-generated catch block
						Alert alert = new Alert(AlertType.ERROR);
						alert.setTitle("Unable to connect to Condor server");
						alert.setHeaderText("Unable to establish connection");
						alert.setContentText(
								"Please check with the helpdesk about the current status of the Condor servers.  If this problem persists, please submit a ticket on this project's repository.");
						e1.printStackTrace();
					}
				}
			}
		});
	}

	public void runBindings() {
		start_run.disableProperty().bind(Bindings
				.when(
						// Requires a run name
						run_name.textProperty()
								// There is some text
								.isNotEqualTo("")
								.and(
										/*
										 * 45 char names max Just makes the
										 * table easier to handle.
										 */
										textCountLimit(run_name, 45))
								.and(textAlphanumeric(run_name)).and(noExistRun(run_name)).and(noExistScript(run_name))
								.and(
										/*
										 * One of the boxes must be check to do
										 * something
										 */
										alignment_check.selectedProperty()
												.or(data_prep_check.selectedProperty()
														.or(counting_check.selectedProperty())
														.or(collect_check.selectedProperty()))))
				.then(false).otherwise(true));
	}

	private String checkSlash(String file_name) {
		if (file_name.length() >= 1) {
			if (file_name.substring(file_name.length() - 1).equals("/")) {
				return file_name;
			} else {
				file_name = file_name + "/";
				return file_name;
			}
		} else {
			return "";
		}
	}

	private BooleanBinding textCountLimit(TextArea run_name, int char_count) {
		BooleanBinding binding = Bindings.createBooleanBinding(() -> (run_name.getText().length() <= char_count),
				run_name.textProperty());
		return binding;
	}

	private BooleanBinding textAlphanumeric(TextArea run_name) {
		BooleanBinding binding = Bindings.createBooleanBinding(
				() -> (!run_name.getText().matches("^.*[^a-zA-Z0-9_].*$")), run_name.textProperty());
		return binding;
	}

	private BooleanBinding noExistRun(TextArea run_name) {
		// Update on run name change
		BooleanBinding bindingRun = Bindings.createBooleanBinding(
				() -> (!new File(checkSlash(destinationDirectory.getText()) + run_name.getText()).isDirectory()),
				run_name.textProperty());
		return bindingRun;
	}

	private BooleanBinding noExistScript(TextArea run_name) {
		// Update on script location change
		BooleanBinding bindingScript = Bindings.createBooleanBinding(
				() -> (!new File(scriptDirectory.getText() + run_name.getText() + ".RData").exists()),
				scriptDirectory.textProperty());
		return bindingScript;
	}

	private RunInstantiator constructRun() {
		// Update no longer in use.
		boolean data_prep = data_prep_check.isSelected();
		boolean alignment = alignment_check.isSelected();
		boolean counting = counting_check.isSelected();
		boolean collect = collect_check.isSelected();
		return new RunInstantiator(data_prep, alignment, counting, collect, run_name.getText());
	}

	private void toggleCountingDep(ToggleGroup group) {
		RSEM.disableProperty().bind(
				Bindings.when(Bowtie.selectedProperty().or(Bowtie2.selectedProperty())).then(false).otherwise(true));
		group.selectedToggleProperty().addListener(new ChangeListener<Toggle>() {
			@Override
			public void changed(ObservableValue<? extends Toggle> observable, Toggle oldValue, Toggle newValue) {
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
		current_tab.selectedIndexProperty().addListener(new ChangeListener<Number>() {
			@Override
			public void changed(ObservableValue<? extends Number> arg0, Number arg1, Number arg2) {
				double currentPercent = ((current_tab.getSelectedIndex() + 2) * percent_per_index) / 100;
				footer_progress.setProgress(currentPercent);
			}
		});
	}

	private void logoutGlow() {
		request.logout();
		Main.setStageTitleDefault();
		Alert alert = new Alert(AlertType.INFORMATION);
		alert.setTitle("Logged out");
		alert.setHeaderText("Logged out");
		alert.setHeaderText("You have logged out of your GLOW session.");

		alert.showAndWait();
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
					if (!"".equals(request.getUsername())) {
						logoutGlow();
						createLoginDialog();
					} else {
						createLoginDialog();
					}
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
		ExtensionFilter R = new ExtensionFilter("R Attribute Files (*.R)", "*.R");
		open_att.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				FileChooser fileChooser = new FileChooser();
				/*
				 * Only looks for R files
				 */
				fileChooser.getExtensionFilters().add(R);
				fileChooser.setSelectedExtensionFilter(R);
				fileChooser.setTitle("Open attribute file");
				File attributeFile = Main.openFileChooser(fileChooser);
				if (attributeFile != null) {
					AttributeAndConfigFileHandler configLoader = new AttributeAndConfigFileHandler();
					try {
						configLoader.setAttributesAttributeFile(attributeFile);
					} catch (FileNotFoundException e) {
						System.out.println("Attribute File not found!");
					}
					String name = attributeFile.getName().split(".R")[0];
					run_name.setText(name);
				}
			}
		});

		template_att.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent event) {
				FileChooser fileChooser = new FileChooser();
				fileChooser.getExtensionFilters().add(R);
				fileChooser.setSelectedExtensionFilter(R);
				fileChooser.setTitle("Open attribute file");
				File att_file = Main.openFileChooser(fileChooser);
				try {
					AttributeAndConfigFileHandler configLoader = new AttributeAndConfigFileHandler();
					configLoader.setAttributesAttributeFile(att_file);
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

				/*
				 * Update absolutely everything.
				 */
				UpdateAttributeSingleton.getInstance().updateAllAttributes(alignment, counting, load, liblistData);
				DirectoryChooser dc = new DirectoryChooser();
				dc.setTitle("Save attribute file");
				File directory = Main.openDirectoryChooser(dc);
				if (directory != null) {
					AttributeFileWriter action = new AttributeFileWriter();
					try {
						action.writeAttributesFile(null, directory.getAbsolutePath());
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
		final String GLOW_USERNAME_PROMPT = "GLOW Username:";
		final String GLOW_PASSWORD_PROMPT = "GLOW Password:";

		// Creates a new dialog which is just a popup
		Dialog<Pair<String, String>> login = new Dialog<>();
		ButtonType loginButtonType = new ButtonType("Login", ButtonData.OK_DONE);
		login.getDialogPane().getButtonTypes().addAll(loginButtonType, ButtonType.CANCEL);
		login.setTitle("Login Dialog");
		login.setHeaderText("GLOW Login");

		// Layout from: http://code.makery.ch/blog/javafx-dialogs-official/
		GridPane grid = new GridPane();
		grid.setHgap(10);
		grid.setVgap(10);
		grid.setPadding(new Insets(20, 150, 10, 10));

		TextField username = new TextField();
		username.setPromptText(GLOW_USERNAME_PROMPT);
		PasswordField password = new PasswordField();
		password.setPromptText(GLOW_PASSWORD_PROMPT);

		grid.add(new Label(GLOW_USERNAME_PROMPT), 0, 0);
		grid.add(username, 1, 0);
		grid.add(new Label(GLOW_PASSWORD_PROMPT), 0, 1);
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
		Button login_button = (Button) login.getDialogPane().lookupButton(loginButtonType);

		login_button.addEventFilter(EventType.ROOT, e -> {
			if (e.getEventType().equals(ActionEvent.ACTION)) {
				e.consume();
				GlowRequest req = new GlowRequest(username.getText(), password.getText());
				request = req;
				login.setTitle("Connecting to Glow Server...Please wait");

				boolean loggedIn = false;
				// Only request cookies when you really need it on windows
				if (runningOnLinux) {
					if (req.requestCookie()) {
						// Login worked, assign it as the current login
						// session
						request = req;
						Task<Object> t = request;
						// Pings the server to verify connection
						new Thread(t).start();
						loggedIn = true;
					}
				} else {
					loggedIn = req.establishSsh() != null;
				}

				if (loggedIn) {
					/*
					 * Give a visual clue at the top of the UI that they are
					 * currently logged into the server. They are actually never
					 * really "logged in", we just keep using their initial
					 * credentials to send further requests and automatically
					 * renew the credentials as necessary.
					 */
					String loggedInHeaderString = "|   GLSeq2 User Interface   -  Currently logged in as "
							+ username.getText() + "   |";
					login.setTitle("Login Dialog");
					Main.setStageTitle(loggedInHeaderString);
					Alert alert = new Alert(AlertType.INFORMATION);
					alert.setTitle("Logged in");
					alert.setHeaderText("You've been logged into GLOW!");
					alert.setContentText(
							"You will be connected to the GLOW server through this client until you either log off explicitly or close the application.");

					alert.showAndWait();

					/*
					 * If the user has logged in and not entered a destination
					 * directory, we'll just auto-populate this as being the
					 * base of their own directory
					 */
					data_request_glow.setText(REQUEST_GLOW);
					ref_request_glow.setText(REQUEST_GLOW);
					if (destinationDirectory.getText() == "") {
						destinationDirectory.setText("/home/GLBRCORG/" + username.getText());
					}

					// Close the dialog as it's job is now complete
					login.close();

					// If the login failed
				} else {
					// http://stackoverflow.com/questions/29911552/how-to-shake-a-login-dialog-in-javafx-8
					try {
						/*
						 * Do a little shake animation. Also clears both the
						 * fields and adds a red border around them to further
						 * indicate that their credentials were wrong.
						 */
						login.setTitle("Login Dialog");
						ShakeTransition anim = new ShakeTransition(login.getDialogPane(), null);
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
		});
	}
}