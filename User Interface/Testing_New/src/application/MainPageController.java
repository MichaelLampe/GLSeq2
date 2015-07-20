package application;

import java.net.URL;
import java.util.ResourceBundle;

import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.RadioButton;
import javafx.scene.control.SingleSelectionModel;
import javafx.scene.control.Slider;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TextArea;
import javafx.scene.control.ToggleGroup;

public class MainPageController implements Initializable {

  @FXML
  Slider gtf_column_slider;
  @FXML
  Label gtf_column_label;

  @FXML
  Slider max_fragment_length_slider;
  @FXML
  Label max_fragment_length_label;

  @FXML
  Slider cores_slider;
  @FXML
  Label cores_label;

  @FXML
  Slider max_buffer_slider;
  @FXML
  Label max_buffer_label;

  @FXML
  Slider dataprep_streams_slider;
  @FXML
  Label dataprep_streams_label;

  @FXML
  Slider alignment_streams_slider;
  @FXML
  Label alignment_streams_label;

  @FXML
  Slider trim_head_slider;
  @FXML
  Label trim_head_label;

  @FXML
  Slider min_read_length_slider;
  @FXML
  Label min_read_length_label;

  @FXML
  ComboBox<String> trim;
  @FXML
  ComboBox<String> alignment_algo;
  @FXML
  ComboBox<String> strandedness;
  @FXML
  ComboBox<String> sequencing_platform;
  @FXML
  ComboBox<String> quality_score;
  @FXML
  ComboBox<String> zipped;
  @FXML
  ComboBox<String> endedness;
  @FXML
  ComboBox<String> presplit;

  @FXML
  TextArea source_directory;
  @FXML
  TextArea destination_directory;
  @FXML
  TextArea countable_sam_directory;
  @FXML
  TextArea previous_run_directory;
  @FXML
  TextArea previous_run_name;
  @FXML
  TextArea strain;
  @FXML
  TextArea reference_genome;
  @FXML
  TextArea reference_fasta;
  @FXML
  TextArea reference_gtf;
  @FXML
  TextArea feature_id_gtf;
  @FXML
  TextArea artificial_fasta;
  @FXML
  TextArea script_directory;
  @FXML
  TextArea trimmomatic_path;
  @FXML
  TextArea fastqc_path;
  @FXML
  TextArea picardTools_directory;
  @FXML
  TextArea bwa_path;
  @FXML
  TextArea bam2wig_path;
  @FXML
  TextArea rockhopper_path;
  @FXML
  TextArea cushaw_gpu_path;
  @FXML
  TextArea tophat_path;
  @FXML
  TextArea cushaw_path;
  @FXML
  TextArea cushaw_index_path;
  @FXML
  TextArea logs_path;
  @FXML
  TextArea run_name;

  @FXML
  CheckBox data_prep_check;
  @FXML
  CheckBox alignment_check;
  @FXML
  CheckBox counting_check;
  @FXML
  CheckBox collect_check;

  @FXML
  Button start_run;

  @FXML
  RadioButton bwa_radio;
  @FXML
  RadioButton bowtie_radio;
  @FXML
  RadioButton bowtie2_radio;
  @FXML
  RadioButton cushaw_radio;
  @FXML
  RadioButton cushaw_gpu_radio;
  @FXML
  RadioButton tophat_radio;
  @FXML
  RadioButton rockhopper_radio;

  @FXML
  TabPane tab_check;
  @FXML
  ProgressBar footer_progress;

  @Override
  public void initialize(URL arg0, ResourceBundle arg1) {
    setupComboBoxes();
    setupStartRun();
    addScrollbarListeners();
    radioButtonGroup();
    progressBarListener();
  }

  private void addScrollbarListeners() {
    // GTF Column
    gtf_column_slider.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int columns = (int) gtf_column_slider.getValue();
        String message;
        if (columns > 1) {
          message = "Your GTF file has " + columns + " columns";
        } else {
          message = "Your GTF file has " + columns + " column";
        }
        gtf_column_label.textProperty().setValue(message);
      }
    });
    // Max Fragment Length
    max_fragment_length_slider.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int max = (int) max_fragment_length_slider.getValue();
        String message;
        message = "The maximal fragment length is " + max;
        max_fragment_length_label.textProperty().setValue(message);
      }
    });
    // Number of Cores
    cores_slider.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int cores = (int) cores_slider.getValue();
        String message;
        if (cores > 1) {
          message = "You are now using " + cores + " cores";
        } else {
          message = "You are now using " + cores + " core";
        }
        cores_label.textProperty().setValue(message);
      }
    });
    // Max CI buffer
    max_buffer_slider.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int buffer = (int) max_buffer_slider.getValue();
        String message;
        message = "Maximum auxiliary buffer when computing CIs is " + buffer;
        max_buffer_label.textProperty().setValue(message);
      }
    });
    // Dataprep Streams
    dataprep_streams_slider.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int streams = (int) dataprep_streams_slider.getValue();
        String message;
        if (streams > 1) {
          message = streams + " computational streams during data preparation";
        } else {
          message = streams + " computational stream during data preparation";
        }
        dataprep_streams_label.textProperty().setValue(message);
      }
    });
    // Alignment Streams
    alignment_streams_slider.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int streams = (int) alignment_streams_slider.getValue();
        String message;
        if (streams > 1) {
          message = streams + " computational streams during alignment";
        } else {
          message = streams + " computational stream during alignment";
        }
        alignment_streams_label.textProperty().setValue(message);
      }
    });
    // Trim Head
    trim_head_slider.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int trim = (int) trim_head_slider.getValue();
        String message;
        if (trim > 1) {
          message = "We are currently trimming " + trim + " base from the head";
        } else {
          message = "We are currently trimming " + trim + " bases from the head";
        }
        trim_head_label.textProperty().setValue(message);
      }
    });
    // Min Read Length
    min_read_length_slider.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int read = (int) min_read_length_slider.getValue();
        String message;
        if (read > 1) {
          message = "Reads under " + read + " bases will be discarded";
        } else {
          message = "Reads under " + read + " base will be discarded";
        }
        min_read_length_label.textProperty().setValue(message);
      }
    });
  }

  private void setupComboBoxes() {
    // Trim Reads
    trim.getItems().addAll("Trim Reads", "Don't Trim Reads");
    trim.setValue("Don't Trim Reads");
    // Strandedness
    strandedness.getItems().addAll("R", "F", "NULL");
    strandedness.setValue("R");
    // Sequencing Platforms
    sequencing_platform.getItems().addAll("illumina", "capillary", "ls454", "solid", "helicos",
        "iontorrent", "pacbio");
    sequencing_platform.setValue("illumina");
    // Quality Score
    quality_score.getItems().addAll("phred33", "phred64");
    quality_score.setValue("phred33");
    // Zipped
    zipped.getItems().addAll("Unzipped", "Zipped");

    // Endedness
    endedness.getItems().addAll("Single", "Paired");

    // Presplit
    presplit.getItems().addAll("Presplit", "Unsplit");
  }

  private void setupStartRun() {
    start_run.setOnAction(new EventHandler<ActionEvent>() {
      @Override
      public void handle(ActionEvent arg0) {
        System.out.println("Start run");
      }
    });
  }

  private void radioButtonGroup() {
    final ToggleGroup group = new ToggleGroup();
    bwa_radio.setToggleGroup(group);
    bowtie_radio.setToggleGroup(group);
    bowtie2_radio.setToggleGroup(group);
    cushaw_radio.setToggleGroup(group);
    cushaw_gpu_radio.setToggleGroup(group);
    tophat_radio.setToggleGroup(group);
    rockhopper_radio.setToggleGroup(group);
  }

  private void progressBarListener() {
    SingleSelectionModel<Tab> current_tab = tab_check.getSelectionModel();
    double percent_per_index = (100.0/tab_check.getTabs().size());
    System.out.println(percent_per_index);
    current_tab.selectedIndexProperty().addListener(new ChangeListener<Number>() {
      @Override
      public void changed(ObservableValue<? extends Number> arg0, Number arg1, Number arg2) {
        double currentPercent = ((current_tab.getSelectedIndex() + 1)*percent_per_index)/100;
        footer_progress.setProgress(currentPercent);
        System.out.println(currentPercent);
      }
    });
  }
}