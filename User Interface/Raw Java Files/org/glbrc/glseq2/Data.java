package org.glbrc.glseq2;

import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.DefaultComboBoxModel;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.JTextPane;
import javax.swing.SpinnerNumberModel;
import javax.swing.UIDefaults;

public class Data extends JDialog {

  private static final long serialVersionUID = 1L;

  private final JTextPane txtchCurrentDataAnd = new JTextPane();
  private final JTextArea txtrawDirectory = new JTextArea(Application.att.getDirectory());
  private final JTextPane txtcDirectoryContainingRaw = new JTextPane();
  private final JPanel panelPreProcessedFiles = new JPanel();
  private final JButton btnPreProcFiles = new JButton("");
  private final JTextArea txtPreProcessedFiles = new JTextArea(Application.att.getDirectoryFq());
  private final JTextPane txtcDirectoryContainingPreprocessed = new JTextPane();
  private final JPanel panelDestinationDir = new JPanel();
  private final JButton btnBaseDir = new JButton("");
  private final JTextArea txtDestinationDirectory = new JTextArea(
      Application.att.getDestinationDirectory());
  private final JTextPane txtcBaseOfDestination = new JTextPane();
  private final JPanel panelRawFiles = new JPanel();
  private final JButton btnRawFiles = new JButton("");
  private final JTextArea txtRawFileNames = new JTextArea(Application.att.getRawFileNames());
  private final JTextPane txtcRawFileNames = new JTextPane();
  private final JButton btnZipped = new JButton();
  private final JButton btnEnded = new JButton();
  private final JPanel panelStrain = new JPanel();
  private final JTextArea txtStrain = new JTextArea(Application.att.getStrain());
  private final JTextPane txtcStrain = new JTextPane();
  private final JPanel panelSubsetOfLibraries = new JPanel();
  private final JTextArea txtSubsetOfLibraries = new JTextArea(Application.att.getLibList());
  private final JTextPane txtcSubsetOfLibraries = new JTextPane();
  private final JPanel panelQualityScores = new JPanel();
  private final JPanel panelSequencing = new JPanel();
  private final JPanel panelStrandedness = new JPanel();
  private final JPanel panelLibraryIdLength = new JPanel();
  private final JTextPane txtcQualityScores = new JTextPane();
  private final JTextPane txtcSequencingPlatforms = new JTextPane();
  private final JTextPane txtcLibraryIdLength = new JTextPane();
  private final JTextPane txtcStrandedness = new JTextPane();
  private final JComboBox<String> comboQualityScores = new JComboBox<String>();
  private final JComboBox<String> comboSequencingPlatforms = new JComboBox<String>();
  private final JComboBox<String> comboStrandedness = new JComboBox<String>();
  private final JSpinner spinLibraryIdLen = new JSpinner();
  private final JPanel panelCountableSam = new JPanel();
  private final JButton btnCountableSam = new JButton("");
  private final JTextArea txtCountableSamDirectory = new JTextArea(
      Application.att.getCountableSamDir());
  private final JTextPane txtcCountableSamFile = new JTextPane();
  private final JButton btnPresplit = new JButton();
  private final JPanel panel = new JPanel();
  private final JButton btnPrevRunDirectory = new JButton("");
  private final JTextArea txtPreviousRunDirectory = new JTextArea(
      Application.att.getPrevRunDirectory());
  private final JTextPane txtcPreviousRunDirectory = new JTextPane();
  private final JPanel panel1 = new JPanel();
  private final JButton btnPrevRunName = new JButton("");
  private final JTextArea txtPreviousRunName = new JTextArea(Application.att.getPrevRunName());
  private final JTextPane txtcPreviousRunName = new JTextPane();

  /**
   * Launch the application.
   */
  public static void main(String[] args) {
    try {
      Data dialog = new Data();
      dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
      dialog.setVisible(true);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * Create the dialog.
   */
  public Data() {
    initGui();
    // Initialize some variables
    if (!Application.att.getQScores().equals("")) {
      comboQualityScores.setSelectedItem(Application.att.getQScores());
    }
    if (!Application.att.getSeqPlatform().equals("")) {
      comboSequencingPlatforms.setSelectedItem(Application.att.getSeqPlatform());
    }
    if (!Application.att.getStrandExtract().equals("")) {
      comboStrandedness.setSelectedItem(Application.att.getLibstrand());
    }
    //
    //
    if (Application.att.getLibNchar().equals("0")) {
      spinLibraryIdLen.setValue(4);
    } else {
      spinLibraryIdLen.setValue(Integer.valueOf(Application.att.getLibNchar()));
    }
    //
    //
    if (Application.att.getUnzipped().equals("TRUE")) {
      btnZipped.setText(ButtonEnums.OptionButton.UNZIPPED.value);
    } else {
      btnZipped.setText(ButtonEnums.OptionButton.ZIPPED.value);
    }
    if (Application.att.getPairedEnd().equals("FALSE")) {
      btnEnded.setText(ButtonEnums.OptionButton.SINGLE.value);
    } else {
      btnEnded.setText(ButtonEnums.OptionButton.PAIRED.value);
    }
    if (Application.att.getPresplit().equals("FALSE")) {
      btnPresplit.setText(ButtonEnums.OptionButton.NO_PRESPLIT.value);
    } else {
      btnPresplit.setText(ButtonEnums.OptionButton.PRESPLIT.value);
    }
  }

  private void initGui() {
    getContentPane().setBackground(Color.LIGHT_GRAY);
    setResizable(false);
    setBounds(100, 100, 1065, 648);
    getContentPane().setLayout(null);
    {
      JPanel buttonPane = new JPanel();
      buttonPane.setBackground(Color.GRAY);
      buttonPane.setBounds(0, 587, 1059, 33);
      buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
      getContentPane().add(buttonPane);
      {
        final JButton okButton = new JButton("Apply and Close");
        okButton.setActionCommand("OK");
        buttonPane.add(okButton);
        okButton.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent arg0) {
            /*
             * Assigns all the values to the attribute file's fields
             */
            Application.att.setDirectory(txtrawDirectory.getText());
            Application.att.setDirectoryFq(txtPreProcessedFiles.getText());
            Application.att.setDestinationDirectory(txtDestinationDirectory.getText());
            Application.att.setRawFileNames(txtRawFileNames.getText());
            Application.att.setCountableSamDir(txtCountableSamDirectory.getText());
            //
            if (btnEnded.getText().equals(ButtonEnums.OptionButton.PAIRED.value)) {
              Application.att.setPairedEnd("TRUE");
            } else {
              Application.att.setPairedEnd("FALSE");
            }
            //
            if (btnZipped.getText().equals(ButtonEnums.OptionButton.ZIPPED.value)) {
              Application.att.setUnzipped("FALSE");
            } else {
              Application.att.setUnzipped("TRUE");
            }
            //
            if (btnPresplit.getText().equals(ButtonEnums.OptionButton.PRESPLIT.value)) {
              Application.att.setPresplit("TRUE");
            } else {
              Application.att.setPresplit("FALSE");
            }
            //
            Application.att.setStrain(txtStrain.getText());
            Application.att.setLibList(txtSubsetOfLibraries.getText());
            Application.att.setQScores(String.valueOf(comboQualityScores.getSelectedItem()));
            Application.att.setSeqPlatform(String.valueOf(comboSequencingPlatforms
                .getSelectedItem()));
            Application.att.setLibstrand(String.valueOf(comboStrandedness.getSelectedItem()));
            Application.att.setLibNchar(String.valueOf(spinLibraryIdLen.getValue()));
            //
            Application.att.setPrevRunDirectory(txtPreviousRunDirectory.getText());
            Application.att.setPrevRunName(txtPreviousRunName.getText());
            //
            dispose();
          }

        });
        getRootPane().setDefaultButton(okButton);
      }
      {
        JButton cancelButton = new JButton("Cancel");
        cancelButton.setActionCommand("Cancel");
        buttonPane.add(cancelButton);
        cancelButton.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent arg0) {
            dispose();
          }
        });
      }
    }
    {
      JPanel panelZipAndStrain = new JPanel();
      panelZipAndStrain.setBackground(Color.LIGHT_GRAY);
      panelZipAndStrain.setBounds(0, 0, 1037, 587);
      getContentPane().add(panelZipAndStrain);
      panelZipAndStrain.setLayout(null);
      {
        txtchCurrentDataAnd.setForeground(Color.DARK_GRAY);
        nimbusFix(Color.LIGHT_GRAY, txtchCurrentDataAnd);
        txtchCurrentDataAnd.setFont(Application.HEADER_FONT);
        txtchCurrentDataAnd.setEditable(false);
        txtchCurrentDataAnd.setText("Current Data and Library Options");
        txtchCurrentDataAnd.setBounds(10, 5, 593, 48);
        panelZipAndStrain.add(txtchCurrentDataAnd);
      }
      {
        JPanel panelRawDirectory = new JPanel();
        panelRawDirectory.setForeground(Color.DARK_GRAY);
        panelRawDirectory.setBackground(Color.LIGHT_GRAY);
        panelRawDirectory.setBounds(10, 60, 653, 42);
        panelZipAndStrain.add(panelRawDirectory);
        panelRawDirectory.setLayout(null);
        {
          JButton btnRawDir = new JButton("");
          btnRawDir.setIcon(new ImageIcon(Data.class
              .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
          btnRawDir.setBounds(271, 0, 38, 37);
          btnRawDir.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
              final JFileChooser chooser = new JFileChooser();
              chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
              chooser.showOpenDialog(txtrawDirectory);
              try {
                File file = chooser.getSelectedFile();
                txtrawDirectory.setText(file.getAbsolutePath());
              } catch (NullPointerException e) {
                // This is meant to be here.
                // It allows for the user to close without selection
              }
            }
          });
          panelRawDirectory.add(btnRawDir);
        }
        txtrawDirectory.setBounds(319, 0, 334, 37);

        panelRawDirectory.add(txtrawDirectory);
        txtcDirectoryContainingRaw.setForeground(Color.DARK_GRAY);
        nimbusFix(Color.LIGHT_GRAY, txtcDirectoryContainingRaw);
        txtcDirectoryContainingRaw.setEditable(false);
        txtcDirectoryContainingRaw.setFont(Application.TEXT_FONT);
        txtcDirectoryContainingRaw.setText("Directory Containing Raw Files");
        txtcDirectoryContainingRaw.setBounds(0, 0, 261, 37);

        panelRawDirectory.add(txtcDirectoryContainingRaw);
      }
      panelPreProcessedFiles.setForeground(Color.DARK_GRAY);
      panelPreProcessedFiles.setBackground(Color.LIGHT_GRAY);
      panelPreProcessedFiles.setLayout(null);
      panelPreProcessedFiles.setBounds(10, 109, 653, 42);

      panelZipAndStrain.add(panelPreProcessedFiles);
      panelDestinationDir.setForeground(Color.DARK_GRAY);
      panelDestinationDir.setBackground(Color.LIGHT_GRAY);
      panelDestinationDir.setLayout(null);
      panelDestinationDir.setBounds(10, 162, 653, 42);

      panelZipAndStrain.add(panelDestinationDir);
      panelRawFiles.setForeground(Color.DARK_GRAY);
      panelRawFiles.setBackground(Color.LIGHT_GRAY);
      panelRawFiles.setLayout(null);
      panelRawFiles.setBounds(10, 215, 653, 42);

      panelZipAndStrain.add(panelRawFiles);
      btnZipped.setFont(new Font("Arial", Font.PLAIN, 15));
      btnZipped.setBounds(690, 342, 337, 59);
      btnZipped.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent arg0) {
          if (btnZipped.getText().equals(ButtonEnums.OptionButton.ZIPPED.value)) {
            btnZipped.setText(ButtonEnums.OptionButton.UNZIPPED.value);
          } else {
            btnZipped.setText(ButtonEnums.OptionButton.ZIPPED.value);
          }
        }
      });

      panelZipAndStrain.add(btnZipped);
      btnEnded.setFont(new Font("Arial", Font.PLAIN, 15));
      btnEnded.setBounds(690, 412, 337, 59);
      btnEnded.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent arg0) {
          if (btnEnded.getText().equals(ButtonEnums.OptionButton.PAIRED.value)) {
            btnEnded.setText(ButtonEnums.OptionButton.SINGLE.value);
          } else {
            btnEnded.setText(ButtonEnums.OptionButton.PAIRED.value);
          }
        }
      });

      panelZipAndStrain.add(btnEnded);
      panelStrain.setForeground(Color.DARK_GRAY);
      panelStrain.setBackground(Color.LIGHT_GRAY);
      panelStrain.setLayout(null);
      panelStrain.setBounds(10, 427, 653, 53);

      panelZipAndStrain.add(panelStrain);
      panelSubsetOfLibraries.setForeground(Color.DARK_GRAY);
      panelSubsetOfLibraries.setBackground(Color.LIGHT_GRAY);
      panelSubsetOfLibraries.setLayout(null);
      panelSubsetOfLibraries.setBounds(10, 491, 653, 53);

      panelZipAndStrain.add(panelSubsetOfLibraries);
      panelQualityScores.setForeground(Color.DARK_GRAY);
      panelQualityScores.setBackground(Color.LIGHT_GRAY);
      panelQualityScores.setBounds(673, 11, 354, 48);

      panelZipAndStrain.add(panelQualityScores);
      panelQualityScores.setLayout(null);
      txtcQualityScores.setForeground(Color.DARK_GRAY);
      nimbusFix(Color.LIGHT_GRAY, txtcQualityScores);
      txtcQualityScores.setEditable(false);
      txtcQualityScores.setFont(Application.TEXT_FONT);
      txtcQualityScores.setText("Quality Scores Format");
      txtcQualityScores.setBounds(10, 0, 334, 20);

      panelQualityScores.add(txtcQualityScores);
      comboQualityScores.setModel(new DefaultComboBoxModel<String>(new String[] { "phred33",
          "phred64" }));
      comboQualityScores.setBounds(10, 28, 334, 20);

      panelQualityScores.add(comboQualityScores);
      panelSequencing.setForeground(Color.DARK_GRAY);
      panelSequencing.setBackground(Color.LIGHT_GRAY);
      panelSequencing.setBounds(673, 123, 354, 48);

      panelZipAndStrain.add(panelSequencing);
      panelSequencing.setLayout(null);
      txtcSequencingPlatforms.setForeground(Color.DARK_GRAY);
      nimbusFix(Color.LIGHT_GRAY, txtcSequencingPlatforms);
      txtcSequencingPlatforms.setEditable(false);
      txtcSequencingPlatforms.setText("Supported Sequencing Platforms");
      txtcSequencingPlatforms.setFont(Application.TEXT_FONT);
      txtcSequencingPlatforms.setBounds(10, 0, 334, 20);

      panelSequencing.add(txtcSequencingPlatforms);
      comboSequencingPlatforms.setModel(new DefaultComboBoxModel<String>(new String[] { "illumina",
          "capillary", "ls454", "solid", "helicos", "iontorrent", "pacbio" }));
      comboSequencingPlatforms.setBounds(10, 28, 334, 20);

      panelSequencing.add(comboSequencingPlatforms);
      panelStrandedness.setForeground(Color.DARK_GRAY);
      panelStrandedness.setBackground(Color.LIGHT_GRAY);
      panelStrandedness.setBounds(673, 64, 354, 48);

      panelZipAndStrain.add(panelStrandedness);
      panelStrandedness.setLayout(null);
      txtcStrandedness.setForeground(Color.DARK_GRAY);
      nimbusFix(Color.LIGHT_GRAY, txtcStrandedness);
      txtcStrandedness.setEditable(false);
      txtcStrandedness.setText("Strandedness of the Library");
      txtcStrandedness.setFont(Application.TEXT_FONT);
      txtcStrandedness.setBounds(10, 0, 334, 20);

      panelStrandedness.add(txtcStrandedness);
      comboStrandedness
          .setModel(new DefaultComboBoxModel<String>(new String[] { "R", "F", "NULL" }));
      comboStrandedness.setBounds(10, 28, 334, 20);

      panelStrandedness.add(comboStrandedness);
      panelLibraryIdLength.setForeground(Color.DARK_GRAY);
      panelLibraryIdLength.setBackground(Color.LIGHT_GRAY);
      panelLibraryIdLength.setBounds(673, 182, 344, 48);

      panelZipAndStrain.add(panelLibraryIdLength);
      panelCountableSam.setLayout(null);
      panelCountableSam.setForeground(Color.DARK_GRAY);
      panelCountableSam.setBackground(Color.LIGHT_GRAY);
      panelCountableSam.setBounds(10, 268, 653, 42);

      panelZipAndStrain.add(panelCountableSam);
      btnPresplit.setFont(new Font("Arial", Font.PLAIN, 15));
      btnPresplit.setBounds(690, 485, 337, 59);
      btnPresplit.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent arg0) {
          if (btnPresplit.getText().equals(ButtonEnums.OptionButton.PRESPLIT.value)) {
            btnPresplit.setText(ButtonEnums.OptionButton.NO_PRESPLIT.value);
          } else {
            btnPresplit.setText(ButtonEnums.OptionButton.PRESPLIT.value);
          }
        }
      });
      panelZipAndStrain.add(btnPresplit);
      panel.setLayout(null);
      panel.setForeground(Color.DARK_GRAY);
      panel.setBackground(Color.LIGHT_GRAY);
      panel.setBounds(10, 321, 653, 42);

      panelZipAndStrain.add(panel);
      panel1.setLayout(null);
      panel1.setForeground(Color.DARK_GRAY);
      panel1.setBackground(Color.LIGHT_GRAY);
      panel1.setBounds(10, 374, 653, 42);

      panelZipAndStrain.add(panel1);
    }
    btnPrevRunName.setIcon(new ImageIcon(Data.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnPrevRunName.setBounds(271, 0, 38, 37);
    btnPrevRunName.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtPreviousRunName);
        try {
          File file = chooser.getSelectedFile();
          txtPreviousRunName.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here
        }
      }
    });

    panel1.add(btnPrevRunName);
    txtPreviousRunName.setBounds(319, 0, 334, 37);

    panel1.add(txtPreviousRunName);
    txtcPreviousRunName.setText("Previous Run Name");
    txtcPreviousRunName.setForeground(Color.DARK_GRAY);
    txtcPreviousRunName.setFont(new Font("Monospaced", Font.PLAIN, 11));
    txtcPreviousRunName.setEditable(false);
    txtcPreviousRunName.setBounds(0, 0, 261, 37);
    nimbusFix(Color.LIGHT_GRAY, txtcPreviousRunName);
    panel1.add(txtcPreviousRunName);
    btnPrevRunDirectory.setIcon(new ImageIcon(Data.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnPrevRunDirectory.setBounds(271, 0, 38, 37);
    btnPrevRunDirectory.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtPreviousRunDirectory);
        try {
          File file = chooser.getSelectedFile();
          txtPreviousRunDirectory.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here
        }
      }
    });

    panel.add(btnPrevRunDirectory);
    txtPreviousRunDirectory.setBounds(319, 0, 334, 37);

    panel.add(txtPreviousRunDirectory);
    txtcPreviousRunDirectory.setText("Previous Run Directory");
    txtcPreviousRunDirectory.setForeground(Color.DARK_GRAY);
    txtcPreviousRunDirectory.setFont(new Font("Monospaced", Font.PLAIN, 11));
    txtcPreviousRunDirectory.setEditable(false);
    txtcPreviousRunDirectory.setBounds(0, 0, 261, 37);
    nimbusFix(Color.LIGHT_GRAY, txtcPreviousRunDirectory);
    panel.add(txtcPreviousRunDirectory);
    btnCountableSam.setIcon(new ImageIcon(Data.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnCountableSam.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent buttionAction) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtCountableSamDirectory);
        try {
          File file = chooser.getSelectedFile();
          txtCountableSamDirectory.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here
        }
      }
    });
    btnCountableSam.setBounds(271, 0, 38, 37);

    panelCountableSam.add(btnCountableSam);
    txtCountableSamDirectory.setBounds(319, 0, 334, 37);

    panelCountableSam.add(txtCountableSamDirectory);
    txtcCountableSamFile.setText("Countable SAM File Directory");
    txtcCountableSamFile.setForeground(Color.DARK_GRAY);
    nimbusFix(Color.LIGHT_GRAY, txtcCountableSamFile);
    txtcCountableSamFile.setFont(Application.TEXT_FONT);
    txtcCountableSamFile.setEditable(false);
    txtcCountableSamFile.setBounds(0, 0, 261, 37);

    panelCountableSam.add(txtcCountableSamFile);
    panelLibraryIdLength.setLayout(null);
    txtcLibraryIdLength.setForeground(Color.DARK_GRAY);
    nimbusFix(Color.LIGHT_GRAY, txtcLibraryIdLength);
    txtcLibraryIdLength.setEditable(false);
    txtcLibraryIdLength.setText("Library ID Length");
    txtcLibraryIdLength.setFont(Application.TEXT_FONT);
    txtcLibraryIdLength.setBounds(0, 11, 141, 37);

    panelLibraryIdLength.add(txtcLibraryIdLength);
    spinLibraryIdLen.setModel(new SpinnerNumberModel(new Integer(4), new Integer(0), null,
        new Integer(1)));
    spinLibraryIdLen.setBounds(169, 0, 61, 48);

    panelLibraryIdLength.add(spinLibraryIdLen);
    txtSubsetOfLibraries.setBounds(319, 11, 334, 31);

    panelSubsetOfLibraries.add(txtSubsetOfLibraries);
    txtcSubsetOfLibraries.setForeground(Color.DARK_GRAY);
    nimbusFix(Color.LIGHT_GRAY, txtcSubsetOfLibraries);
    txtcSubsetOfLibraries.setEditable(false);
    txtcSubsetOfLibraries.setText("Subset of Libraries to Process");
    txtcSubsetOfLibraries.setFont(Application.TEXT_FONT);
    txtcSubsetOfLibraries.setBounds(0, 11, 309, 31);

    panelSubsetOfLibraries.add(txtcSubsetOfLibraries);
    txtStrain.setBounds(319, 11, 334, 31);

    panelStrain.add(txtStrain);
    txtcStrain.setForeground(Color.DARK_GRAY);
    nimbusFix(Color.LIGHT_GRAY, txtcStrain);
    txtcStrain.setEditable(false);
    txtcStrain.setText("Strain");
    txtcStrain.setFont(Application.TEXT_FONT);
    txtcStrain.setBounds(0, 11, 309, 31);

    panelStrain.add(txtcStrain);
    btnRawFiles.setIcon(new ImageIcon(Data.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnRawFiles.setBounds(271, 0, 38, 37);
    btnRawFiles.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtRawFileNames);
        try {
          File file = chooser.getSelectedFile();
          txtRawFileNames.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here
        }
      }
    });

    panelRawFiles.add(btnRawFiles);
    txtRawFileNames.setBounds(319, 0, 334, 37);
    panelRawFiles.add(txtRawFileNames);
    txtcRawFileNames.setForeground(Color.DARK_GRAY);
    nimbusFix(Color.LIGHT_GRAY, txtcRawFileNames);
    txtcRawFileNames.setEditable(false);
    txtcRawFileNames.setText("Raw File Names");
    txtcRawFileNames.setFont(Application.TEXT_FONT);
    txtcRawFileNames.setBounds(0, 0, 261, 37);

    panelRawFiles.add(txtcRawFileNames);
    btnBaseDir.setIcon(new ImageIcon(Data.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnBaseDir.setBounds(271, 0, 38, 37);
    btnBaseDir.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtDestinationDirectory);
        try {
          File file = chooser.getSelectedFile();
          txtDestinationDirectory.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here
        }
      }
    });

    panelDestinationDir.add(btnBaseDir);
    txtDestinationDirectory.setBounds(319, 0, 334, 37);

    panelDestinationDir.add(txtDestinationDirectory);
    txtcBaseOfDestination.setForeground(Color.DARK_GRAY);
    nimbusFix(Color.LIGHT_GRAY, txtcBaseOfDestination);
    txtcBaseOfDestination.setEditable(false);
    txtcBaseOfDestination.setText("Base of Destination Directory");
    txtcBaseOfDestination.setFont(Application.TEXT_FONT);
    txtcBaseOfDestination.setBounds(0, 0, 261, 37);

    panelDestinationDir.add(txtcBaseOfDestination);
    btnPreProcFiles.setIcon(new ImageIcon(Data.class
        .getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
    btnPreProcFiles.setBounds(271, 0, 38, 37);
    btnPreProcFiles.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        final JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.showOpenDialog(txtPreProcessedFiles);
        try {
          File file = chooser.getSelectedFile();
          txtPreProcessedFiles.setText(file.getAbsolutePath());
        } catch (NullPointerException e) {
          // This is meant to be here
        }
      }
    });

    panelPreProcessedFiles.add(btnPreProcFiles);
    txtPreProcessedFiles.setBounds(319, 0, 334, 37);

    panelPreProcessedFiles.add(txtPreProcessedFiles);
    txtcDirectoryContainingPreprocessed.setForeground(Color.DARK_GRAY);
    nimbusFix(Color.LIGHT_GRAY, txtcDirectoryContainingPreprocessed);
    txtcDirectoryContainingPreprocessed.setEditable(false);
    txtcDirectoryContainingPreprocessed.setText("Directory Containing Pre-Processed Files");
    txtcDirectoryContainingPreprocessed.setFont(new Font("Courier New", Font.PLAIN, 11));
    txtcDirectoryContainingPreprocessed.setBounds(0, 0, 261, 37);

    panelPreProcessedFiles.add(txtcDirectoryContainingPreprocessed);
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
