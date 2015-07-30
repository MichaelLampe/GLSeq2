package application;

import java.io.File;
import java.util.List;

import javafx.concurrent.Task;

public class RunWorker extends Task<Object> {
  private List<String> args;
  private RunTableEntry currentRun;

  public RunWorker(List<String> args) {
    this.args = args;
    this.currentRun = null;
  }

  public RunWorker(List<String> args, RunTableEntry currentRun) {
    this.args = args;
    this.currentRun = currentRun;
  }

  @Override
  protected Object call() throws Exception {
    // Script generated based on arguments from the RunOptions class.
    // It calls the R script with the correct user parameters
    // Build process w/ args again
    String scriptDirectory = Attributes.getInstance().attributesCollection.get("base.dir")
        .getValue();
    try {
      Process process = new ProcessBuilder(args).directory(new File(scriptDirectory)).inheritIO().start();
      // Wait until it is done or it crashes.
      process.waitFor();
      if (process.exitValue() != 0){
        currentRun.setError();
      } else{
        currentRun.setComplete();
      }
    } catch (Exception e) {
      if (currentRun != null) {
        currentRun.setError();
      }
    }
    return null;
  }
}
