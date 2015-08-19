package application;

import java.io.IOException;
import java.net.URL;
import java.util.List;
import java.util.ResourceBundle;

import javafx.beans.binding.Bindings;
import javafx.beans.binding.BooleanBinding;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.concurrent.Task;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.fxml.Initializable;
import javafx.scene.control.CheckBox;
import javafx.scene.control.SingleSelectionModel;
import javafx.scene.control.Tab;
import javafx.scene.control.TextArea;
import javafx.scene.control.ToggleGroup;
import javafx.scene.control.cell.PropertyValueFactory;

public final class MainPageController extends MainPageItems implements Initializable {

  public static ToggleGroup group;

  @Override
  public void initialize(URL arg0, ResourceBundle arg1) {
    setupComboBoxes();
    setupStartRun();
    addScrollbarListeners();
    group = radioButtonGroup();
    progressBarListener();
    runBindings();
  }

  /*
   * All scroll bars have listeners that allow for dynamic text update to their
   * associated labels. Those are all created here.
   */
  private void addScrollbarListeners() {
    // GTF Column
    gtfFeatureColumn.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int columns = (int) gtfFeatureColumn.getValue();
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
    fragMaxLength.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int max = (int) fragMaxLength.getValue();
        String message;
        message = "The maximal fragment length is " + max;
        max_fragment_length_label.textProperty().setValue(message);
      }
    });
    // Number of Cores
    numberCores.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int cores = (int) numberCores.getValue();
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
    ciMem.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int buffer = (int) ciMem.getValue();
        String message;
        message = "Maximum auxiliary buffer when computing CIs is " + buffer;
        max_buffer_label.textProperty().setValue(message);
      }
    });
    // Dataprep Streams
    numberStreamsDataPrep.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int streams = (int) numberStreamsDataPrep.getValue();
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
    numberStreams.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int streams = (int) numberStreams.getValue();
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
    trimHead.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int trim = (int) trimHead.getValue();
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
    minTrim.valueProperty().addListener(new ChangeListener<Object>() {
      @Override
      public void changed(ObservableValue<?> arg0, Object arg1, Object arg2) {
        int read = (int) minTrim.getValue();
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
    seqPlatform.getItems().addAll("illumina", "capillary", "ls454", "solid", "helicos",
        "iontorrent", "pacbio");
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
        AttributeActions action = new AttributeActions();
        UpdateAttribute.getInstance().updateAttributes();
        String attributeFileLocation = null;
        try {
          if (run_name.getText().equals("")) {
            attributeFileLocation = action.writeAttributesFile(null, null);
          } else {
            attributeFileLocation = action.writeAttributesFile(run_name.getText(), null);
          }
          action.saveConfigFile("AttributesConfig.txt");
        } catch (IOException e) {
          e.printStackTrace();
        }
        // Creates a run instance
        Run currentRun = constructRun();
        // Grabs the args from that instance
        // If true, run on HTcondor, else do not
        List<String> args = currentRun.constructArgs(htcondor_check.isSelected(),
            attributeFileLocation);
        if (args != null) {
          // Might add some logic in the future that starts an SSH on
          // windows/mac
          Task<Object> task = new RunWorker(args, currentRun.getActiveRun());
          // Run on a different thread so it doesn't lock up the UI.
          new Thread(task).start();
        } else {
          System.out.println("Please select something for the script to do.");
        }
      }
    });
  }

  public void runBindings() {
    start_run.disableProperty().bind(
        Bindings
            .when(
            // Requires a run name
                run_name
                    .textProperty()
                    // There is some text
                    .isNotEqualTo("")
                    .and(
                    // 45 char names max
                    // Just makes the table easier to handle.
                        textCountLimit(run_name, 45))
                    .and(
                    // One of the boxes must be check to do something.
                        alignment_check.selectedProperty().or(
                            data_prep_check.selectedProperty()
                                .or(counting_check.selectedProperty())
                                .or(collect_check.selectedProperty())))).then(false)
            .otherwise(true));
  }

  private BooleanBinding textCountLimit(TextArea run_name, int char_count) {
    BooleanBinding binding = Bindings.createBooleanBinding(
        () -> (run_name.getText().length() <= char_count), run_name.textProperty());
    return binding;
  }

  private Run constructRun() {
    // Update no longer in use.
    boolean data_prep = data_prep_check.isSelected();
    boolean alignment = alignment_check.isSelected();
    boolean counting = counting_check.isSelected();
    boolean collect = collect_check.isSelected();
    return new Run(data_prep, alignment, counting, collect, run_name.getText());
  }

  /*
   * Alignment radio buttons need to be grouped into toggle groups to allow only
   * one selection. They are put into one of those here.
   */
  private ToggleGroup radioButtonGroup() {
    final ToggleGroup group = new ToggleGroup();
    BWA.setToggleGroup(group);
    Bowtie.setToggleGroup(group);
    Bowtie2.setToggleGroup(group);
    Cushaw.setToggleGroup(group);
    Cushaw_Gpu.setToggleGroup(group);
    TopHat.setToggleGroup(group);
    Rockhopper.setToggleGroup(group);
    BWA.setSelected(true);
    return group;
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
        double currentPercent = ((current_tab.getSelectedIndex() + 1) * percent_per_index) / 100;
        footer_progress.setProgress(currentPercent);
      }
    });
  }
}