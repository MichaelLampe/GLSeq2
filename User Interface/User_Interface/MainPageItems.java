package application;

import javafx.fxml.FXML;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.control.MenuItem;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.RadioButton;
import javafx.scene.control.Slider;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.TextArea;
import javafx.scene.control.Tooltip;
import javafx.scene.control.TreeView;

class MainPageItems {
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
  protected Slider trimHead;
  @FXML
  protected Label trim_head_label;

  @FXML
  protected Slider minTrim;
  @FXML
  protected Label min_read_length_label;

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
  protected TextArea strain;
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
  protected CheckBox RSEM;
  @FXML
  protected CheckBox HTseq;
  @FXML
  protected CheckBox FeatureCounts;
  @FXML
  protected CheckBox Cufflinks;
  @FXML
  protected CheckBox rockhopperCount;

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
  protected TabPane tab_check;
  @FXML
  protected ProgressBar footer_progress;

  @FXML
  protected TableView<RunTableEntry> table_runs;
  @FXML
  protected TableColumn<RunTableEntry, String> column_run_name;
  @FXML
  protected TableColumn<RunTableEntry, String> column_status;
  @FXML
  protected TableColumn<RunTableEntry, String> column_start;
  @FXML
  protected TableColumn<RunTableEntry, String> column_end;
  @FXML
  protected TableColumn<RunTableEntry, String> column_duration;

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
}
