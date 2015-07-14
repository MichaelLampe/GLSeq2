package org.glbrc.glseq2.Project_Files;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPasswordField;
import javax.swing.JTextField;

import com.jcabi.ssh.SSHByPassword;
import com.jcabi.ssh.Shell;

public class ShellSsh extends JDialog {
  // This class name is funny.
  private static final long serialVersionUID = 1L;

  private String address;
  private int port;
  private String login;
  private String password;
  private String[] serverCommands;
  private final JTextField txtServer = new JTextField();
  private final JLabel txtcServer = new JLabel(
      "Server (Host Name or IP) ('GLBRC' will utilize Scarcity Servers)");
  private final JLabel txtcPort = new JLabel("Port Number (Default of 22)");
  private final JTextField txtPort = new JTextField();
  private final JLabel txtcUsername = new JLabel("Username");
  private final JTextField txtUsername = new JTextField();
  private final JLabel txtcPassword = new JLabel("Password");
  private final JPasswordField txtPassword = new JPasswordField();
  private final JLabel lblNewLabel_1 = new JLabel("");
  private final JButton btnRunScript = new JButton("Run Scripts");

  /**
   * Constructs a panel that allows login credential input
   * 
   * @param address
   *          host name of IP address of a server
   * @param port
   *          The port to use
   * @param login
   *          The login of the user
   * @param password
   *          The password of the user
   */
  public ShellSsh() {
    initGui();
  }

  /**
   * GUI.
   */
  public void initGui() {
    setTitle("GLBRC SSH Login");
    
    setResizable(false);
    setBounds(100, 100, 442, 293);
    
    
    txtServer.setText("GLBRC");
    txtServer.setBounds(10, 36, 414, 20);
    txtServer.setColumns(10);
    getContentPane().setLayout(null);

    getContentPane().add(txtServer);
    txtcServer.setBounds(10, 11, 414, 20);

    getContentPane().add(txtcServer);
    txtcPort.setBounds(10, 61, 414, 20);

    getContentPane().add(txtcPort);
    txtPort.setText("22");
    txtPort.setColumns(10);
    txtPort.setBounds(10, 86, 414, 20);

    getContentPane().add(txtPort);
    txtcUsername.setBounds(10, 111, 414, 20);

    getContentPane().add(txtcUsername);
    txtUsername.setColumns(10);
    txtUsername.setBounds(10, 136, 414, 20);

    getContentPane().add(txtUsername);
    txtcPassword.setBounds(10, 161, 414, 20);

    getContentPane().add(txtcPassword);
    txtPassword.setColumns(10);
    txtPassword.setBounds(10, 186, 414, 20);

    getContentPane().add(txtPassword);
    lblNewLabel_1.setBounds(10, 11, 46, 14);

    getContentPane().add(lblNewLabel_1);
    btnRunScript.setBounds(10, 217, 414, 34);

    getContentPane().add(btnRunScript);

    /*
     * Action Listeners
     */
    btnRunScript.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        updateFields();
        submitJobs();
      }
    });
  }

  private void updateFields() {
    address = txtServer.getText();
    port = Integer.valueOf(txtPort.getText());
    login = txtUsername.getText();
    password = String.valueOf(txtPassword.getPassword());
  }

  /**
   * The top level executive function that forces the shell to execute.
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
      shell = new SSHByPassword(server, port, login, password);
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
   * @param commands
   *          - An array containing all of the shell commands
   * @param servers
   *          - A list of all of the servers to use.
   */
  private void runStack(String[] commands, ArrayList<String> servers) {
    for (int i = 0; i < commands.length; i++) {
      try {
        new Shell.Plain(new SSHByPassword(servers.get(i), port, login, password)).exec(commands[i]);
      } catch (UnknownHostException e) {
        System.out.println("Unknown host.  Please check servers.");
      } catch (IOException e) {
        System.out.println("Input output problem, please check the commands that were sent.");
      }
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
      Shell shell = new SSHByPassword(server, port, login, password);
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