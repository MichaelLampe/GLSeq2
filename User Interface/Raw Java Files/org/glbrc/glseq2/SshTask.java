package org.glbrc.glseq2;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.SwingWorker;

import com.jcabi.ssh.SSHByPassword;
import com.jcabi.ssh.Shell;

// Extends Swing worker to control GUI lockup upon running (Aka it leaves the GUI all fine).
public class SshTask extends SwingWorker<List<Integer>, Integer> {

  private final String address;
  private final int port;
  private final String username;
  private final String password;
  private String[] serverCommands;

  /**
   * Constructor for a single SSH task. Creates a new login event and executes
   * the task
   * 
   * @param address
   *          - Server address
   * @param port
   *          - Port number
   * @param username
   *          - Personal user name
   * @param password
   *          - Adjacent Password
   */
  public SshTask(String address, int port, String username, String password) {
    this.address = address;
    this.port = port;
    this.username = username;
    this.password = password;
  }

  @Override
  protected List<Integer> doInBackground() throws Exception {
    submitJobs();
    return null;
  }

  /**
   * Server ping function that allows outside sources to access the availability
   * of the server for this login information.
   */
  public boolean ping() {
    String temp;
    if (address.equals("GLBRC")) {
      temp = "scarcity-1.glbrc.org";
    } else {
      temp = address;
    }
    try {
      new Shell.Plain(new SSHByPassword(temp, port, username, password))
          .exec("echo 'ping from GLBRC_UI'");

    } catch (UnknownHostException e) {
      return false;
    } catch (IOException e) {
      return false;
    }
    return true;
  }

  /**
   * The top level executive function that forces the shell to execute.
   * 
   * <p>
   * Handles some special shortcut cases for our servers to allow batch
   * scheduling (Servers are treated as a list as opposed to a single string).
   */
  private void submitJobs() {
    ArrayList<String> servers = new ArrayList<String>();
    if (address.equals("GLBRC")) {
      servers = glbrcServers();
    } else {
      servers.add(address);
    }
    // Each index here represents a single server that will have a single line
    // of commands executed.
    Application.updating("Constructing server commands");
    serverCommands = new String[servers.size()];
    int server = 0;
    for (QueuedRun q : Application.tabsRun.getQueues()) {
      String scriptDir = q.getSelectedRun().getScriptDirectory(
          q.getSelectedAttributes().getScriptDirectory());
      String move = "cd " + scriptDir;
      if (serverCommands[server] == null) {
        // Changes directories to the script one where the RScript will be run.
        serverCommands[server] = move;
      } else {
        serverCommands[server] += " ; " + move;
      }
      // Special case, GPU accel must run on our server as server 10 which has a
      // GPU
      if (q.getSelectedAttributes().getaAlgor().equals("Cushaw")) {
        if (q.getSelectedAttributes().getGpuAccel().equals("TRUE")) {
          // The list is made counting backwards so that in order we will use
          // the fastest computers (And 10 has a GPU, so that's index 0)
          String attFilePath = startAnalysis(q, servers.get(0));
          q.getSelectedRun().setAttributeFilePath(attFilePath);
          List<String> args = q.getSelectedRun().returnArgs();
          String command = " && ";
          for (String arg : args) {
            command += arg + " ";
          }
          serverCommands[0] += command;
          // Only move to next server if we are at server 10 currently.
          if (server == 0) {
            server++;
          }
        }
      } else {
        String attFilePath = startAnalysis(q, servers.get(server));
        q.getSelectedRun().setAttributeFilePath(attFilePath);
        List<String> args = q.getSelectedRun().returnArgs();
        String command = " && ";
        for (String arg : args) {
          command += arg + " ";
        }
        serverCommands[server] += command;
        server++;
      }
      if (server >= serverCommands.length) {
        server = 0;
      }
    }
    Application.updating("Running commands via SSH.");
    // Remove all the tabs and then run it
    Application.tabsRun.removeQueues();
    runStack(serverCommands, servers);
  }

  /**
   * Generates all the GLBRC server names.
   * 
   * @return Severs Returns a list of all the GLBRC servers.
   */
  private ArrayList<String> glbrcServers() {
    // Set our servers as default
    ArrayList<String> servers = new ArrayList<String>();
    for (int i = 10; i > 0; i--) {
      servers.add("scarcity-" + i + ".glbrc.org");
    }
    return servers;
  }

  /**
   * This script transfers an attribute file to a server via SSH, creates a
   * directory for the file, and then executes the Rscript.
   * 
   * @param QueuedRun
   *          A single queued run
   */
  private String startAnalysis(QueuedRun que, String server) {
    String dir = setupEnvironment(server);
    String fileName = que.getSelectedRun().getAttributeFilePath();
    String attributeFileLocation = dir + "/" + que.attributeFileName();
    // Linux command cat > {file} which writes a file to the server when the
    // File is used as the input stream.
    String script = "cat > " + attributeFileLocation;
    copyAttributeFiles(fileName, script, server);
    return attributeFileLocation;
  }

  /**
   * Copies the attribute file from the native Windows platform onto the Linux
   * link.
   * 
   * @param file
   *          The file name that will be transfered
   * @param script
   *          The command that will be executed to
   * @param server
   *          Which server to connect to
   * 
   */
  private void copyAttributeFiles(String file, String script, String server) {
    Shell shell = null;
    try {
      // Access a server via SSH w/ login credentials
      shell = new SSHByPassword(server, port, username, password);
    } catch (UnknownHostException e1) {
      System.out.println("Unable to connect to host.");
      return;
    }

    File attFile = new File(file);
    try {
      new Shell.Safe(shell).exec(script, new FileInputStream(attFile), null, null);
    } catch (FileNotFoundException e) {
      System.out.println("Could not find attribute file to transfer to server.");
    } catch (IOException e) {
      System.out.println("Error when undergoing transfer and connecting to server. "
          + " Please check all your credentials and the server.");
    } catch (NullPointerException e) {
      // This is called when shell for some reason is null (Aka it failed. That
      // should be never as the catch returns out of the method.
    }
  }

  /**
   * For all the commands, where each command represents one server, launch the
   * shell and run it on the server.
   * 
   * <p>
   * This uses some parallelization, but due to the need for Java to hold onto
   * the thread Larger runs where the number of jobs exceeds the cores on the
   * submitting computer Will notice that the maximum number of concurrent jobs
   * is equal to the number of cores.
   * 
   * <p>
   * Parallel for loop Modified from from
   * http://stackoverflow.com/a/4010275/4978569
   * 
   * @param commands
   *          - An array containing all of the shell commands
   * @param servers
   *          - A list of all of the servers to use.
   */
  private void runStack(String[] commands, ArrayList<String> servers) {
    // Sets the threads that can be used as all available threads.
    ExecutorService exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
    int at = 0;
    try {
      for (final String command : commands) {
        // Keeps track of the server. Must be final.
        final int cur = at;
        exec.submit(new Runnable() {
          @Override
          public void run() {
            try {
              // Don't bother sending commands that have no use (Length of 0)
              if (command.length() > 0) {
                System.out.println("On server " + cur + " we are running the following command:");
                System.out.println(command);

                // new Shell.Plain(new SSHByPassword(servers.get(cur), port,
                // username, password))
                // .exec(command);
              }
            } catch (Exception e) {
              System.out.println("ok");
            }
          }
        });
        // Make sure to iterate the to keep the server track going!
        at++;
      }
    } finally {
      exec.shutdown();
    }
  }

  /**
   * Makes a standard folder in the user's starting path that stores their
   * attribute files.
   * 
   * @param server
   *          Which server to connect to via SSH
   */
  private String setupEnvironment(String server) {
    String directory = null;
    try {
      Shell shell = new SSHByPassword(server, port, username, password);
      directory = new Shell.Plain(shell).exec("pwd");
      directory = directory.replace("\n", "");
      directory = directory + "/AttributeFiles";
      // Tries to create the directory for the user.
      new Shell.Plain(shell).exec("mkdir " + directory);
    } catch (IOException e) {
      e.printStackTrace();
      System.out.println("Error executing the command or connecting to the host.");
    }
    return directory;
  }
}
