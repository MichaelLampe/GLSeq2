package application;

import java.io.File;
import java.io.IOException;
import java.util.List;

import javafx.concurrent.Task;

public class RunWorker extends Task<Object> {
  List<String> args;

  public RunWorker(List<String> args) {
    this.args = args;
  }

  @Override
  protected Object call() {
    // Script generated based on arguments from the RunOptions class.
    // It calls the R script with the correct user parameters
    // Build process w/ args again
    String scriptDirectory = Attributes.getInstance().attributesCollection.get("base.dir")
        .getValue();
    try {
      // Builds a new process
      // Changes directory to the script dir
      // Changes the IO stream to the current java instance (System's)
      // Finally starts the process
      Process process = new ProcessBuilder(args).directory(new File(scriptDirectory)).inheritIO().start();
      try {
        process.waitFor();
        System.out.println("Script exited with value:" + process.exitValue());
        // Make sure that the process is dead.
        process.destroyForcibly();
      } catch (InterruptedException e) {
        System.out.println("Script was interrupted during run.");
      }
    } catch (IOException e1) {
      e1.printStackTrace();
      System.out.println("Problem starting script");
    }
    return null;
  }
}
