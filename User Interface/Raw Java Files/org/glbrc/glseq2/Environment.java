package org.glbrc.glseq2;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextPane;
import javax.swing.UIDefaults;
import javax.swing.border.EmptyBorder;

public class Environment extends JDialog {

  private static final long serialVersionUID = 1L;
  private final JPanel contentPanel = new JPanel();

  private final JTextPane txtchEnvironmentOptions = new JTextPane();
  /*
   * / GLSEQ Script
   */
  private final JTextArea txtGlSeqDirectory = new JTextArea();
  private final JTextPane txtcGlSeqDirectory = new JTextPane();
  private final JButton btnGlSeqDirectory = new JButton("");
  /*
   * Picard tools menu
   */
  private final JTextArea txtPicard = new JTextArea();
  private final JPanel panelPicard = new JPanel();
  private final JTextPane txtcPicard = new JTextPane();
  private final JButton btnPicard = new JButton("");
  /*
   * Fastqc menu
   */
  private final JTextArea txtFastqc = new JTextArea();
  private final JPanel panelFastqc = new JPanel();
  private final JTextPane txtcFastqc = new JTextPane();
  private final JButton btnFastqc = new JButton("");
  /*
   * Bam2wig menu
   */
  private final JTextArea txtBam2Wig = new JTextArea();
  private final JPanel panelBam2Wig = new JPanel();
  private final JTextPane txtcBam2Wig = new JTextPane();
  private final JButton btnBam2Wig = new JButton("");
  /*
   * Bwa menu
   */
  private final JPanel panelBwa = new JPanel();
  private final JTextPane txtcBwa = new JTextPane();
  private final JTextArea txtBwa = new JTextArea();
  private final JButton btnBwa = new JButton("");
  /*
   * Trimmomatic menu
   */
  private final JPanel panelTrimmomatic = new JPanel();
  private final JTextArea txtTrimmomatic = new JTextArea();
  private final JTextPane txtcTrimmomatic = new JTextPane();
  private final JButton btnTrimmomatic = new JButton("");
  /*
   * Cushaw menu
   */
  private final JTextArea txtCushaw = new JTextArea();
  private final JPanel panelCushaw = new JPanel();
  private final JTextPane txtcCushaw = new JTextPane();
  /*
   * Cushaw index menu
   */
  private final JPanel panelCushawIndex = new JPanel();
  private final JTextPane txtcCushawIndex = new JTextPane();
  private final JTextArea txtCushawIndex = new JTextArea();
  private final JButton btnCushawIndex = new JButton("");
  private final JButton btnCushaw = new JButton("");
  /*
   * Cushaw GPU menu
   */
  private final JPanel panelCushawGpu = new JPanel();
  private final JTextPane txtcCushawGpu = new JTextPane();
  private final JTextArea txtCushawGpu = new JTextArea();
  private final JButton btnCushawGpu = new JButton("");
  /*
   * Movement options
   */
  private final JPanel panelButton = new JPanel();
  private final JButton okButton = new JButton("Apply and Close");
  private final JButton cancelButton = new JButton("Cancel");
  private final JPanel panel = new JPanel();
  private final JScrollPane scrollPane = new JScrollPane();
  private final JScrollPane scrollPane1 = new JScrollPane();
  private final JScrollPane scrollPane2 = new JScrollPane();
  private final JScrollPane scrollPane3 = new JScrollPane();
  private final JScrollPane scrollPane4 = new JScrollPane();
  private final JScrollPane scrollPane5 = new JScrollPane();
  private final JScrollPane scrollPane6 = new JScrollPane();
  private final JScrollPane scrollPane7 = new JScrollPane();
  private final JScrollPane scrollPane8 = new JScrollPane();
  private final JPanel tophatPanel = new JPanel();
  private final JTextPane txtpnPathToTophat = new JTextPane();
  private final JButton btnTopHat = new JButton("");
  private final JScrollPane scrollPane9 = new JScrollPane();
  private final JTextArea txtTopHat = new JTextArea();
  private final JPanel verboseLogPanel = new JPanel();
  private final JTextPane txtcVerboseLog = new JTextPane();
  private final JButton btnVerboseLog = new JButton("");
  private final JScrollPane verboseScroll = new JScrollPane();
  private final JTextArea txtVerboseLog = new JTextArea();

  /**
   * Launch the application.
   */
  public static void main(String[] args) {
    try {
      Environment dialog = new Environment();
      dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
      dialog.setVisible(true);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * Create the dialog.
   */
  public Environment() {
    setResizable(false);
    initGui();
    btnTrimmomatic.setIcon(new ImageIcon(Environment.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnTrimmomatic.setBounds(191, 12, 38, 37);
    btnTrimmomatic.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtTrimmomatic);
        try {
          File file = chooser.getSelectedFile();
          txtTrimmomatic.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here.
        }
      }
    });

    panelTrimmomatic.add(btnTrimmomatic);
    scrollPane1.setBounds(239, 11, 151, 37);

    panelTrimmomatic.add(scrollPane1);
    scrollPane1.setViewportView(txtTrimmomatic);
    // Initialize some variables
    txtTrimmomatic.setText(Application.att.getTrimPath());
    btnPicard.setIcon(new ImageIcon(Environment.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnPicard.setBounds(191, 12, 38, 37);
    btnPicard.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtPicard);
        try {
          File file = chooser.getSelectedFile();
          txtPicard.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here.
        }
      }
    });

    panelPicard.add(btnPicard);
    scrollPane2.setBounds(239, 11, 151, 37);

    panelPicard.add(scrollPane2);
    scrollPane2.setViewportView(txtPicard);
    txtPicard.setText(Application.att.getPicardToolsPath());
    btnFastqc.setIcon(new ImageIcon(Environment.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnFastqc.setBounds(191, 12, 38, 37);
    btnFastqc.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtFastqc);
        try {
          File file = chooser.getSelectedFile();
          txtFastqc.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here.
        }
      }
    });

    panelFastqc.add(btnFastqc);
    scrollPane3.setBounds(239, 11, 151, 37);

    panelFastqc.add(scrollPane3);
    scrollPane3.setViewportView(txtFastqc);
    txtFastqc.setText(Application.att.getFastqcPath());
    btnBwa.setIcon(new ImageIcon(Environment.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnBwa.setBounds(191, 12, 38, 37);
    btnBwa.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtBwa);
        try {
          File file = chooser.getSelectedFile();
          txtBwa.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here.
        }
      }
    });

    panelBwa.add(btnBwa);
    scrollPane4.setBounds(239, 11, 151, 37);

    panelBwa.add(scrollPane4);
    scrollPane4.setViewportView(txtBwa);
    txtBwa.setText(Application.att.getBwaPath());
    btnBam2Wig.setIcon(new ImageIcon(Environment.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnBam2Wig.setBounds(191, 12, 38, 37);
    btnBam2Wig.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtBam2Wig);
        try {
          File file = chooser.getSelectedFile();
          txtBam2Wig.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here.
        }
      }
    });

    panelBam2Wig.add(btnBam2Wig);
    scrollPane5.setBounds(239, 11, 151, 37);

    panelBam2Wig.add(scrollPane5);
    scrollPane5.setViewportView(txtBam2Wig);
    txtBam2Wig.setText(Application.att.getBam2WigPath());
    btnCushaw.setIcon(new ImageIcon(Environment.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnCushaw.setBounds(191, 12, 38, 37);
    btnCushaw.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtCushaw);
        try {
          File file = chooser.getSelectedFile();
          txtCushaw.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here.
        }
      }
    });

    panelCushaw.add(btnCushaw);
    scrollPane7.setBounds(239, 11, 151, 37);

    panelCushaw.add(scrollPane7);
    scrollPane7.setViewportView(txtCushaw);
    txtCushaw.setText(Application.att.getCushawPath());
    btnCushawIndex.setIcon(new ImageIcon(Environment.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnCushawIndex.setBounds(191, 12, 38, 37);
    btnCushawIndex.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtCushawIndex);
        try {
          File file = chooser.getSelectedFile();
          txtCushawIndex.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here.
        }
      }
    });

    panelCushawIndex.add(btnCushawIndex);
    scrollPane6.setBounds(239, 11, 151, 37);

    panelCushawIndex.add(scrollPane6);
    scrollPane6.setViewportView(txtCushawIndex);
    txtCushawIndex.setText(Application.att.getCushawIndexPath());
    btnCushawGpu.setIcon(new ImageIcon(Environment.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnCushawGpu.setBounds(191, 11, 38, 37);
    btnCushawGpu.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtCushawGpu);
        try {
          File file = chooser.getSelectedFile();
          txtCushawGpu.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here.
        }
      }
    });

    panelCushawGpu.add(btnCushawGpu);
    scrollPane8.setBounds(239, 11, 151, 37);

    panelCushawGpu.add(scrollPane8);
    scrollPane8.setViewportView(txtCushawGpu);
    txtCushawGpu.setText(Application.att.getCushawGpuPath());
    panel.setLayout(null);
    panel.setForeground(Color.DARK_GRAY);
    panel.setBackground(Color.LIGHT_GRAY);
    panel.setBounds(10, 48, 400, 59);

    contentPanel.add(panel);
    txtcGlSeqDirectory.setText("Path to GLSeq Scripts");
    txtcGlSeqDirectory.setForeground(Color.DARK_GRAY);
    nimbusFix(Color.LIGHT_GRAY, txtcGlSeqDirectory);
    txtcGlSeqDirectory.setFont(new Font("Monospaced", Font.PLAIN, 11));
    txtcGlSeqDirectory.setEditable(false);
    txtcGlSeqDirectory.setBounds(10, 11, 171, 37);

    panel.add(txtcGlSeqDirectory);
    btnGlSeqDirectory.setIcon(new ImageIcon(Environment.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnGlSeqDirectory.setBounds(191, 11, 38, 37);
    btnGlSeqDirectory.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtGlSeqDirectory);
        try {
          File file = chooser.getSelectedFile();
          txtGlSeqDirectory.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here.
        }
      }
    });
    panel.add(btnGlSeqDirectory);
    scrollPane.setBounds(239, 11, 151, 37);

    panel.add(scrollPane);
    scrollPane.setViewportView(txtGlSeqDirectory);
    txtGlSeqDirectory.setText(Application.att.getScriptDirectory());
    tophatPanel.setLayout(null);
    tophatPanel.setForeground(Color.DARK_GRAY);
    tophatPanel.setBackground(Color.LIGHT_GRAY);
    tophatPanel.setBounds(413, 48, 400, 59);

    contentPanel.add(tophatPanel);
    txtpnPathToTophat.setText("Path to TopHat");
    txtpnPathToTophat.setForeground(Color.DARK_GRAY);
    txtpnPathToTophat.setFont(new Font("Monospaced", Font.PLAIN, 11));
    txtpnPathToTophat.setEditable(false);
    txtpnPathToTophat.setBounds(10, 11, 171, 37);

    tophatPanel.add(txtpnPathToTophat);
    btnTopHat.setIcon(new ImageIcon(Environment.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnTopHat.setBounds(191, 11, 38, 37);
    btnTopHat.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtTopHat);
        try {
          File file = chooser.getSelectedFile();
          txtTopHat.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here.
        }
      }
    });
    tophatPanel.add(btnTopHat);
    scrollPane9.setBounds(239, 11, 151, 37);

    tophatPanel.add(scrollPane9);
    txtTopHat.setText(Application.att.getTopHatPath());

    scrollPane9.setViewportView(txtTopHat);
    verboseLogPanel.setLayout(null);
    verboseLogPanel.setForeground(Color.DARK_GRAY);
    verboseLogPanel.setBackground(Color.LIGHT_GRAY);
    verboseLogPanel.setBounds(413, 526, 400, 59);
    
    contentPanel.add(verboseLogPanel);
    txtcVerboseLog.setText("Verbose Test Log");
    txtcVerboseLog.setForeground(Color.DARK_GRAY);
    txtcVerboseLog.setFont(new Font("Monospaced", Font.PLAIN, 11));
    txtcVerboseLog.setEditable(false);
    txtcVerboseLog.setBounds(10, 11, 171, 37);
    
    verboseLogPanel.add(txtcVerboseLog);
    btnVerboseLog.setIcon(new ImageIcon(Environment.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnVerboseLog.setBounds(191, 11, 38, 37);
    btnVerboseLog.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtTrimmomatic);
        try {
          File file = chooser.getSelectedFile();
          txtVerboseLog.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here.
        }
      }
    });
    verboseLogPanel.add(btnVerboseLog);
    verboseScroll.setBounds(239, 11, 151, 37);
    
    verboseLogPanel.add(verboseScroll);
    nimbusFix(Color.LIGHT_GRAY, txtcVerboseLog);
    txtVerboseLog.setText(Application.att.getDestDirTest());
    verboseScroll.setViewportView(txtVerboseLog);

  }

  private void initGui() {
    setBackground(Color.LIGHT_GRAY);
    setBounds(100, 100, 829, 663);
    getContentPane().setLayout(new BorderLayout());
    contentPanel.setBackground(Color.LIGHT_GRAY);
    contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    getContentPane().add(contentPanel, BorderLayout.CENTER);
    contentPanel.setLayout(null);
    {
      panelFastqc.setForeground(Color.DARK_GRAY);
      panelFastqc.setBackground(Color.LIGHT_GRAY);
      panelFastqc.setLayout(null);
      panelFastqc.setBounds(10, 236, 400, 59);
      contentPanel.add(panelFastqc);
      {
        txtcFastqc.setForeground(Color.DARK_GRAY);
        nimbusFix(Color.LIGHT_GRAY, txtcFastqc);
        txtcFastqc.setText("Path to Fastqc");
        txtcFastqc.setFont(Application.TEXT_FONT);
        txtcFastqc.setEditable(false);
        txtcFastqc.setBounds(10, 11, 171, 37);
        panelFastqc.add(txtcFastqc);
      }
    }
    {
      panelPicard.setForeground(Color.DARK_GRAY);
      panelPicard.setBackground(Color.LIGHT_GRAY);
      panelPicard.setLayout(null);
      panelPicard.setBounds(10, 177, 400, 59);
      contentPanel.add(panelPicard);
      {
        txtcPicard.setForeground(Color.DARK_GRAY);
        nimbusFix(Color.LIGHT_GRAY, txtcPicard);
        txtcPicard.setText("Path to PicardTools Jar Directory");
        txtcPicard.setFont(Application.TEXT_FONT);
        txtcPicard.setEditable(false);
        txtcPicard.setBounds(10, 11, 171, 37);
        panelPicard.add(txtcPicard);
      }
    }
    {
      panelTrimmomatic.setForeground(Color.DARK_GRAY);
      panelTrimmomatic.setBackground(Color.LIGHT_GRAY);
      panelTrimmomatic.setLayout(null);
      panelTrimmomatic.setBounds(10, 118, 400, 59);
      contentPanel.add(panelTrimmomatic);
      {
        txtcTrimmomatic.setForeground(Color.DARK_GRAY);
        nimbusFix(Color.LIGHT_GRAY, txtcTrimmomatic);
        txtcTrimmomatic.setText("Path to Trimmomatic");
        txtcTrimmomatic.setFont(Application.TEXT_FONT);
        txtcTrimmomatic.setEditable(false);
        txtcTrimmomatic.setBounds(10, 11, 171, 37);
        panelTrimmomatic.add(txtcTrimmomatic);
      }
    }
    {
      panelCushaw.setForeground(Color.DARK_GRAY);
      panelCushaw.setBackground(Color.LIGHT_GRAY);
      panelCushaw.setLayout(null);
      panelCushaw.setBounds(10, 467, 400, 59);
      contentPanel.add(panelCushaw);
      {
        txtcCushaw.setForeground(Color.DARK_GRAY);
        nimbusFix(Color.LIGHT_GRAY, txtcCushaw);
        txtcCushaw.setText("Path to CUSHAW");
        txtcCushaw.setFont(Application.TEXT_FONT);
        txtcCushaw.setEditable(false);
        txtcCushaw.setBounds(10, 11, 171, 37);
        panelCushaw.add(txtcCushaw);
      }
    }
    {
      panelCushawGpu.setForeground(Color.DARK_GRAY);
      panelCushawGpu.setBackground(Color.LIGHT_GRAY);
      panelCushawGpu.setLayout(null);
      panelCushawGpu.setBounds(10, 526, 400, 59);
      contentPanel.add(panelCushawGpu);
      {
        txtcCushawGpu.setForeground(Color.DARK_GRAY);
        nimbusFix(Color.LIGHT_GRAY, txtcCushawGpu);
        txtcCushawGpu.setText("Path to CUSHAW-GPU");
        txtcCushawGpu.setFont(Application.TEXT_FONT);
        txtcCushawGpu.setEditable(false);
        txtcCushawGpu.setBounds(10, 11, 171, 37);
        panelCushawGpu.add(txtcCushawGpu);
      }
    }
    {
      panelBam2Wig.setForeground(Color.DARK_GRAY);
      panelBam2Wig.setBackground(Color.LIGHT_GRAY);
      panelBam2Wig.setLayout(null);
      panelBam2Wig.setBounds(10, 352, 400, 59);
      contentPanel.add(panelBam2Wig);
      {
        txtcBam2Wig.setForeground(Color.DARK_GRAY);
        nimbusFix(Color.LIGHT_GRAY, txtcBam2Wig);
        txtcBam2Wig.setText("Path to Bam2Wig Shell Script");
        txtcBam2Wig.setFont(Application.TEXT_FONT);
        txtcBam2Wig.setEditable(false);
        txtcBam2Wig.setBounds(10, 11, 171, 37);
        panelBam2Wig.add(txtcBam2Wig);
      }
    }
    {
      panelBwa.setForeground(Color.DARK_GRAY);
      panelBwa.setBackground(Color.LIGHT_GRAY);
      panelBwa.setLayout(null);
      panelBwa.setBounds(10, 295, 400, 59);
      contentPanel.add(panelBwa);
      {
        txtcBwa.setForeground(Color.DARK_GRAY);
        nimbusFix(Color.LIGHT_GRAY, txtcBwa);
        txtcBwa.setText("Path to BWA");
        txtcBwa.setFont(Application.TEXT_FONT);
        txtcBwa.setEditable(false);
        txtcBwa.setBounds(10, 11, 171, 37);
        panelBwa.add(txtcBwa);
      }
    }
    {
      txtchEnvironmentOptions.setFont(Application.HEADER_FONT);
      txtchEnvironmentOptions.setForeground(Color.DARK_GRAY);
      nimbusFix(Color.LIGHT_GRAY, txtchEnvironmentOptions);
      txtchEnvironmentOptions.setText("Environment Options");
      txtchEnvironmentOptions.setEditable(false);
      txtchEnvironmentOptions.setBounds(10, 4, 496, 33);
      contentPanel.add(txtchEnvironmentOptions);
    }
    {
      panelCushawIndex.setLayout(null);
      panelCushawIndex.setForeground(Color.DARK_GRAY);
      panelCushawIndex.setBackground(Color.LIGHT_GRAY);
      panelCushawIndex.setBounds(10, 408, 400, 59);
      contentPanel.add(panelCushawIndex);
      {
        txtcCushawIndex.setText("Path to CUSHAW Index");
        txtcCushawIndex.setForeground(Color.DARK_GRAY);
        nimbusFix(Color.LIGHT_GRAY, txtcCushawIndex);
        txtcCushawIndex.setFont(Application.TEXT_FONT);
        txtcCushawIndex.setEditable(false);
        txtcCushawIndex.setBounds(10, 11, 171, 37);
        panelCushawIndex.add(txtcCushawIndex);
      }
    }
    {
      panelButton.setBackground(Color.GRAY);
      panelButton.setLayout(new FlowLayout(FlowLayout.RIGHT));
      getContentPane().add(panelButton, BorderLayout.SOUTH);
      {
        okButton.setActionCommand("OK");
        okButton.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent arg0) {
            Application.att.setTrimPath(txtTrimmomatic.getText());
            Application.att.setPicardToolsPath(txtPicard.getText());
            Application.att.setFastqcPath(txtFastqc.getText());
            Application.att.setBwaPath(txtBwa.getText());
            Application.att.setBam2WigPath(txtBam2Wig.getText());
            Application.att.setCushawPath(txtCushaw.getText());
            Application.att.setCushawIndexPath(txtCushawIndex.getText());
            Application.att.setCushawGpuPath(txtCushawGpu.getText());
            Application.att.setScriptDirectory(txtGlSeqDirectory.getText());
            Application.att.setTopHatPath(txtTopHat.getText());
            Application.att.setDestDirTest(txtVerboseLog.getText());
            dispose();
          }
        });
        panelButton.add(okButton);
        getRootPane().setDefaultButton(okButton);
      }
      {
        cancelButton.setActionCommand("Cancel");
        panelButton.add(cancelButton);
        cancelButton.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent arg0) {
            dispose();
          }
        });
      }
    }
  }

  // Fixes a bug in Nimbus where it overrides the desired JTextPane
  // Background colors
  void nimbusFix(Color background, JTextPane pane) {
    UIDefaults defaults = new UIDefaults();
    defaults.put("TextPane[Enabled].backgroundPainter", background);
    pane.putClientProperty("Nimbus.Overrides", defaults);
    pane.putClientProperty("Nimbus.Overrides.InheritDefaults", true);
    pane.setBackground(background);
  }
}
