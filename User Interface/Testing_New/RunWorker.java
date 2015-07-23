package application;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
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
    ProcessBuilder script = new ProcessBuilder(args);
    String scriptDirectory = Attributes.getInstance().attributesCollection.get("base.dir").getValue();
    script.directory(new File(scriptDirectory));
    Process process = null;
    try {
      process = script.start();
      final InputStream is = process.getInputStream();
      System.out.println(is);
      final InputStream es = process.getErrorStream();
      new Thread(new Runnable() {
        public void run() {
          try {
            BufferedReader reader = new BufferedReader(new InputStreamReader(is));
            String line;
            while ((line = reader.readLine()) != null) {
              System.out.println(line);
            }
            reader = new BufferedReader(new InputStreamReader(es));
            while ((line = reader.readLine()) != null) {
              System.out.println(line);
            }
          } catch (IOException e) {
            e.printStackTrace();
          } finally {
            try {
              is.close();
            } catch (IOException e) {
              e.printStackTrace();
            }
          }
        }
      }).start();
      System.out.println("Started script");
      // Script is done when it says done
      try {
        process.waitFor();
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
