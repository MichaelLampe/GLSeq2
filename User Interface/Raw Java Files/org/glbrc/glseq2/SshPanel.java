package org.glbrc.glseq2;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPasswordField;
import javax.swing.JTextField;
import javax.swing.border.EmptyBorder;

public class SshPanel extends JDialog {

  private static final long serialVersionUID = 1L;

  private final JPanel contentPanel = new JPanel();

  private final JLabel txtcServer = new JLabel(
      "Server (Host Name or IP) ('GLBRC' will utilize Scarcity Servers)");
  private final JLabel txtcPort = new JLabel("Port Number (Default of 22)");
  private final JLabel txtcUsername = new JLabel("Username");
  private final JLabel txtcPassword = new JLabel("Password");
  private final JButton btnRunScripts = new JButton("Run Scripts");

  private final JTextField txtServer;
  private final JTextField txtPort;
  private final JTextField txtUsername;
  private final JPasswordField txtPassword;

  /**
   * Launch the application.
   */
  public static void main(String[] args) {
    try {
      SshPanel dialog = new SshPanel();
      dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
      dialog.setVisible(true);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * Create the dialog.
   */
  public SshPanel() {
    txtServer = new JTextField();
    txtPort = new JTextField();
    txtUsername = new JTextField();
    txtPassword = new JPasswordField();
    initGui();
  }

  /**
   * Initializes the GUI for the JDialog frame.
   * 
   */
  public void initGui() {
    setTitle("GLBRC SSH Login");
    setBounds(100, 100, 458, 341);
    getContentPane().setLayout(new BorderLayout());
    contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    getContentPane().add(contentPanel, BorderLayout.CENTER);
    contentPanel.setLayout(null);
    {
      txtcServer.setBounds(10, 10, 414, 20);
      contentPanel.add(txtcServer);
    }
    {
      txtServer.setText("GLBRC");
      txtServer.setColumns(10);
      txtServer.setBounds(10, 35, 414, 30);
      contentPanel.add(txtServer);
    }
    {
      txtcPort.setBounds(10, 70, 414, 20);
      contentPanel.add(txtcPort);
    }
    {
      txtPort.setText("22");
      txtPort.setColumns(10);
      txtPort.setBounds(10, 95, 414, 30);
      contentPanel.add(txtPort);
    }
    {
      txtcUsername.setBounds(10, 130, 426, 20);
      contentPanel.add(txtcUsername);
    }
    {
      txtUsername.setColumns(10);
      txtUsername.setBounds(10, 155, 414, 30);
      contentPanel.add(txtUsername);
    }
    {
      txtcPassword.setBounds(10, 190, 475, 20);
      contentPanel.add(txtcPassword);
    }
    {
      txtPassword.setColumns(10);
      txtPassword.setBounds(10, 215, 414, 30);
      contentPanel.add(txtPassword);
    }
    {
      btnRunScripts.setBounds(10, 256, 414, 34);
      contentPanel.add(btnRunScripts);
    }

    btnRunScripts.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        String server = txtServer.getText();
        int port = Integer.valueOf(txtPort.getText());
        String username = txtUsername.getText();
        String password = String.valueOf(txtPassword.getPassword());
        SshTask ssh = new SshTask(server, port, username, password);
        boolean auth = ssh.ping();
        if (auth) {
          ssh.execute();
          Application.updating("Username and password verified.  Submitting jobs to server(s)");
          setVisible(false);
        } else {
          JOptionPane.showMessageDialog(contentPanel,
              "Error with login. Please check credentials and server address.");
          txtPassword.setText("");
        }
      }
    });
  }
}