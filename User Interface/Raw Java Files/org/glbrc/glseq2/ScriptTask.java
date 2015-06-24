package org.glbrc.glseq2;

///////////////////////////////////////////////////////////////////////////////
//                   
// Main Class File:  GLSeq2_Main_Application.java
// File:             ScriptTask.java
//
//
// Author:           Michael Lampe | MrLampe@Wisc.edu
///////////////////////////////////////////////////////////////////////////////

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import javax.swing.SwingWorker;

// SwingWorker allows for parallel execution
public class ScriptTask extends SwingWorker<List<Integer>, Integer> {
  private Run localRun;
  private Attributes localAttributes;
  private boolean batch;

  /**
   * This constructor is called when one run is done. It has the added bonus of
   * also reporting a few more statistics that would get cumbersome (currently)
   * with ~12 runs going on.
   * 
   * @param localRun
   *          - The run given to it that it will use to generate the args
   * @param localAttributes
   *          - The attributes it will use to save a copy.
   */
  public ScriptTask(Run localRun, Attributes localAttributes) {
    this.localRun = localRun;
    this.localAttributes = localAttributes;
    batch = false;
  }

  /**
   * This constructor utilizes all the tabbed runs to generate a queue of runs.
   * They will run sequentially and non-dependently via the ';' in the Linux
   * command line.
   * 
   */
  public ScriptTask() {
    batch = true;
    // No declaration. This means we'll be doing a batch run here.

  }

  @Override
  protected List<Integer> doInBackground() throws IOException {
    startScript();
    return null;
  }

  private void startScript() {
    /*
     * Saves both the attribute and run config file. This is just pulling data
     * from the correct sources to get a good file path for the saved files.
     */
    if (!batch) {
      try {
        localAttributes
            .saveConfigFile(localRun.getDestinationDirectory(localAttributes
                .getDestinationDirectory())
                + "/AttributesConfigFile_"
                + localAttributes.runId
                + ".txt");
        localRun.saveConfigFile(localRun.getDestinationDirectory(localAttributes
            .getDestinationDirectory()) + "/RunConfigFile_" + localAttributes.runId + ".txt");
      } catch (IOException e) {
        // This is meant to be here
      }
      // Indicate to user that the script has started
    }
    List<Run> allRun = new ArrayList<Run>();
    List<Attributes> attributes = new ArrayList<Attributes>();
    if (batch) {
      List<QueuedRun> allRuns = Application.tabsRun.getQueues();
      //
      for (QueuedRun run : allRuns) {
        allRun.add(run.getSelectedRun());
        attributes.add(run.getSelectedAttributes());
      }
      //
      for (QueuedRun run : allRuns) {
        Application.tabsRun.removeQueue(run);
        QueuedRun.count--;
      }
    } else {
      allRun.add(localRun);
      //
      attributes.add(localAttributes);
    }
    // Script generated based on arguments from the RunOptions class.
    // It calls the R script with the correct user parameters
    for (int i = 0; i < allRun.size(); i++) {
      ProcessBuilder script = new ProcessBuilder(allRun.get(i).returnArgs());
      script.directory(new File(allRun.get(i).getScriptDirectory(
          attributes.get(i).getScriptDirectory())));
      Process process = null;
      try {
        process = script.start();
        System.out.println("Started script");
        final InputStream is = process.getInputStream();
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
        // Error message if the script is incorrectly run. Did not select
        // directory or options are easy reasons for this error.
      } catch (IOException e) {
        Application.updating("Sorry, there was a problem running the script. "
            + " Please check that all the options have been correctly selected.");
        Application.updating(String.valueOf(e));
      }
      try {
        // Script is done when it says done
        process.waitFor();
        Application.updating("The Top Script has now completely executed. "
            + "Processing will continue in your destination folder(s) until complete.");
      } catch (InterruptedException e) {
        Application.updating("Error running the script. Script interrupted.");
      }
    }
  }
}
