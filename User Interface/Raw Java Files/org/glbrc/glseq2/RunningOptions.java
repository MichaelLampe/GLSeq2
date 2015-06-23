package org.glbrc.glseq2;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.JTextPane;
import javax.swing.SpinnerNumberModel;
import javax.swing.UIDefaults;
import javax.swing.border.EmptyBorder;

public class RunningOptions extends JDialog {
  /*
   * Various variables regarding the UI of the Script Running
   */
  private static final long serialVersionUID = 1L;
  private final JPanel contentPanel = new JPanel();
  //
  private final JPanel panelManualHolder = new JPanel();
  //
  private final JTextArea txtProtocolId = new JTextArea();
  private final JPanel panelProtocol = new JPanel();
  private final JTextPane txtcProtocolId = new JTextPane();
  //
  private final JSpinner spinParallelExpression = new JSpinner();
  private final JPanel panelParallelExpressionComp = new JPanel();
  private final JTextPane txtcParallelExpressionComp = new JTextPane();
  //
  private final JPanel panelParallelDataComp = new JPanel();
  private final JTextPane txtParallelDataPrep = new JTextPane();
  private final JSpinner spinParallelDataPrep = new JSpinner();
  //
  private final JTextPane txtcMaximalFragmentLength = new JTextPane();
  private final JPanel panelMaxFragLen = new JPanel();
  private final JSpinner spinMaxFragLen = new JSpinner();
  //
  private final JSpinner spinMaxBuffer = new JSpinner();
  private final JPanel panelMaxBufferConfInts = new JPanel();
  private final JTextPane txtcMaximumAuxiliaryBuffer = new JTextPane();
  //
  private final JButton btnExtractCoverage = new JButton();
  private final JButton btnOutputGenomeBam = new JButton();
  //
  private final JPanel panelAlignment = new JPanel();
  private final JTextPane txtcAlignmentAlgorithm = new JTextPane();
  private final JComboBox<String> comboAlignmentAlgo = new JComboBox<String>();
  //
  private final JPanel panelCounting = new JPanel();
  private final JCheckBox checkFeatureCount = new JCheckBox("Feature Counts");
  private final JCheckBox checkHtSeq = new JCheckBox("HTSeq");
  private final JCheckBox checkRsem = new JCheckBox("RSEM");
  private final JTextPane txtcCountingMethods = new JTextPane();
  //
  private final JTextPane txtchPipelineOptions = new JTextPane();
  //
  private final JPanel panelNumberOfCores = new JPanel();
  private final JTextPane txtcNumberOfCores = new JTextPane();
  private final JSpinner spinNumberOfCores = new JSpinner();
  //
  private final JPanel buttonPane = new JPanel();
  private final JPanel panel = new JPanel();
  private final JPanel panel1 = new JPanel();
  //
  private final JPanel panelReferenceGenome = new JPanel();
  private final JTextArea txtReferenceGenome = new JTextArea();
  private final JTextPane txtcReferenceGenome = new JTextPane();
  //
  private final JPanel panelReferenceFasta = new JPanel();
  private final JTextArea txtReferenceFasta = new JTextArea();
  private final JTextPane txtcReferenceFasta = new JTextPane();
  //
  private final JPanel panelReferenceFeatures = new JPanel();
  private final JTextArea txtReferenceFeatures = new JTextArea();
  private final JTextPane txtcReferenceFeatures = new JTextPane();
  //
  private final JPanel panelReferenceId = new JPanel();
  private final JTextArea txtReferenceId = new JTextArea();
  private final JTextPane txtcReferenceId = new JTextPane();
  //
  private final JTextPane txtcColumnGtf = new JTextPane();
  private final JSpinner spinColumnGtf = new JSpinner();
  //
  private final JButton okButton = new JButton("Apply and Close");
  private final JButton cancelButton = new JButton("Cancel");
  //
  private final JPanel panelGlow = new JPanel();
  private final JTextPane txtpnGlowLoginGoing = new JTextPane();

  /**
   * Launch the application.
   */
  public static void main(String[] args) {
    try {
      RunningOptions dialog = new RunningOptions();
      dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
      dialog.setVisible(true);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * Create the dialog.
   */
  public RunningOptions() {
    initGui();
    //
    spinNumberOfCores.setValue(Integer.valueOf(Application.att.getCores()));
    spinParallelExpression.setValue(Integer.valueOf(Application.att.getStreams()));
    spinParallelDataPrep.setValue(Integer.valueOf(Application.att.getStreamsDataPrep()));
    spinMaxFragLen.setValue(Integer.valueOf(Application.att.getFragMaxLength()));
    spinMaxBuffer.setValue(Integer.valueOf(Application.att.getCiMem()));
    txtReferenceGenome.setText(Application.att.getRGenome());
    txtReferenceFasta.setText(Application.att.getFasta());
    txtReferenceFeatures.setText(Application.att.getGff());
    txtReferenceId.setText(Application.att.getIdAttr());
    spinColumnGtf.setValue(Integer.valueOf(Application.att.getFeatureColumn()));
    //
    if (!String.valueOf(comboAlignmentAlgo.getSelectedItem()).contains("Bowtie")) {
      checkRsem.setSelected(false);
      checkRsem.setEnabled(false);
    } else {
      checkRsem.setEnabled(true);
    }
    //
    if (String.valueOf(comboAlignmentAlgo.getSelectedItem()).contains("Cushaw")) {
      if (Application.att.getGpuAccel().equals("TRUE")) {
        comboAlignmentAlgo.setSelectedItem("Cushaw-GPU");
      }
    }
    //
    if (Application.att.getFeatureCounts().contains("FeatureCounts")) {
      checkFeatureCount.setSelected(true);
    }
    //
    if (Application.att.getHtseq().contains("HTSeq")) {
      checkHtSeq.setSelected(true);
    }
    //
    if (Application.att.getRsem().contains("RSEM")
        && Application.att.getaAlgor().contains("Bowtie")) {
      checkRsem.setSelected(true);
    }
    //
    {
      panel.setLayout(null);
      panel.setForeground(Color.DARK_GRAY);
      panel.setBackground(Color.LIGHT_GRAY);
      panel.setBounds(0, 516, 695, 59);
      panelManualHolder.add(panel);
    }
    {
      panelReferenceGenome.setLayout(null);
      panelReferenceGenome.setForeground(Color.DARK_GRAY);
      panelReferenceGenome.setBackground(Color.LIGHT_GRAY);
      panelReferenceGenome.setBounds(10, 0, 983, 40);
      panelManualHolder.add(panelReferenceGenome);
    }
    {
      txtReferenceGenome.setBounds(254, 0, 719, 37);
      panelReferenceGenome.add(txtReferenceGenome);
    }
    {
      txtcReferenceGenome.setText("Reference Genome");
      txtcReferenceGenome.setForeground(Color.DARK_GRAY);
      txtcReferenceGenome.setFont(new Font("Monospaced", Font.PLAIN, 11));
      txtcReferenceGenome.setEditable(false);
      nimbusFix(Color.LIGHT_GRAY, txtcReferenceGenome);
      txtcReferenceGenome.setBounds(0, 0, 235, 37);
      panelReferenceGenome.add(txtcReferenceGenome);
    }
    {
      panelReferenceFasta.setLayout(null);
      panelReferenceFasta.setBackground(Color.LIGHT_GRAY);
      panelReferenceFasta.setBounds(10, 41, 983, 40);
      panelManualHolder.add(panelReferenceFasta);
    }
    {
      txtReferenceFasta.setBounds(255, 0, 718, 37);
      panelReferenceFasta.add(txtReferenceFasta);
    }
    {
      txtcReferenceFasta.setText("Name of Reference FASTA File");
      txtcReferenceFasta.setForeground(Color.DARK_GRAY);
      txtcReferenceFasta.setFont(new Font("Monospaced", Font.PLAIN, 11));
      nimbusFix(Color.LIGHT_GRAY, txtcReferenceFasta);
      txtcReferenceFasta.setEditable(false);
      txtcReferenceFasta.setBounds(0, 0, 235, 37);
      panelReferenceFasta.add(txtcReferenceFasta);
    }
    {
      panelReferenceFeatures.setLayout(null);
      panelReferenceFeatures.setForeground(Color.DARK_GRAY);
      panelReferenceFeatures.setBackground(Color.LIGHT_GRAY);
      panelReferenceFeatures.setBounds(10, 81, 983, 40);
      panelManualHolder.add(panelReferenceFeatures);
    }
    {
      txtReferenceFeatures.setBounds(255, 0, 718, 37);
      panelReferenceFeatures.add(txtReferenceFeatures);
    }
    {
      txtcReferenceFeatures.setText("Name of Reference Genomic Features File");
      txtcReferenceFeatures.setForeground(Color.DARK_GRAY);
      nimbusFix(Color.LIGHT_GRAY, txtcReferenceFeatures);
      txtcReferenceFeatures.setFont(new Font("Monospaced", Font.PLAIN, 11));
      txtcReferenceFeatures.setEditable(false);
      txtcReferenceFeatures.setBounds(0, 0, 235, 37);
      panelReferenceFeatures.add(txtcReferenceFeatures);
    }
    {
      panelReferenceId.setLayout(null);
      panelReferenceId.setForeground(Color.DARK_GRAY);
      panelReferenceId.setBackground(Color.LIGHT_GRAY);
      panelReferenceId.setBounds(10, 122, 983, 47);
      panelManualHolder.add(panelReferenceId);
    }
    {
      txtReferenceId.setBounds(255, 0, 718, 37);
      panelReferenceId.add(txtReferenceId);
    }
    {
      txtcReferenceId.setText("Reference Feature ID");
      txtcReferenceId.setForeground(Color.DARK_GRAY);
      nimbusFix(Color.LIGHT_GRAY, txtcReferenceId);
      txtcReferenceId.setFont(new Font("Monospaced", Font.PLAIN, 11));
      txtcReferenceId.setEditable(false);
      txtcReferenceId.setBounds(0, 0, 235, 37);
      panelReferenceId.add(txtcReferenceId);
    }
    txtcColumnGtf.setText("Number of Columns in the GTF File");
    txtcColumnGtf.setForeground(Color.DARK_GRAY);
    nimbusFix(Color.LIGHT_GRAY, txtcColumnGtf);
    txtcColumnGtf.setEditable(false);
    txtcColumnGtf.setBounds(526, 384, 368, 40);
    //
    panelManualHolder.add(txtcColumnGtf);
    spinColumnGtf.setBounds(905, 384, 77, 40);
    //
    panelManualHolder.add(spinColumnGtf);
    {
      panel1.setLayout(null);
      panel1.setBackground(Color.LIGHT_GRAY);
      panel1.setBounds(314, 0, 699, 180);
      contentPanel.add(panel1);
    }
    txtProtocolId.setText(Application.run.getProtocolId());
    comboAlignmentAlgo.setSelectedItem(Application.att.getaAlgor());
    //
    if (Application.att.getStrandExtract().equals("TRUE")) {
      btnExtractCoverage.setText(ButtonEnums.OptionButton.EXTRACT.value);
    } else {
      btnExtractCoverage.setText(ButtonEnums.OptionButton.NO_EXTRACT.value);
    }
    //
    if (Application.att.getGenoBam().equals("TRUE")) {
      btnOutputGenomeBam.setText(ButtonEnums.OptionButton.OUTPUT.value);
    } else {
      btnOutputGenomeBam.setText(ButtonEnums.OptionButton.NO_OUTPUT.value);
    }
  }

  private void initGui() {
    setBounds(100, 100, 1029, 693);
    getContentPane().setLayout(new BorderLayout());
    contentPanel.setBackground(Color.LIGHT_GRAY);
    contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    getContentPane().add(contentPanel, BorderLayout.CENTER);
    contentPanel.setLayout(null);
    {
      panelGlow.setBackground(Color.ORANGE);
      panelGlow.setLayout(null);
      panelGlow.setBounds(0, 0, 298, 180);
      contentPanel.add(panelGlow);
    }
    {
      txtpnGlowLoginGoing.setText("Glow Login Going Here");
      txtpnGlowLoginGoing.setBounds(39, 38, 204, 20);
      panelGlow.add(txtpnGlowLoginGoing);
    }
    {
      panelManualHolder.setBackground(Color.LIGHT_GRAY);
      panelManualHolder.setBounds(10, 191, 1003, 525);
      contentPanel.add(panelManualHolder);
      panelManualHolder.setLayout(null);
      {
        panelProtocol.setBounds(295, 11, 394, 59);
        panel1.add(panelProtocol);
        panelProtocol.setForeground(Color.DARK_GRAY);
        panelProtocol.setBackground(Color.LIGHT_GRAY);
        panelProtocol.setLayout(null);
        {
          txtProtocolId.setBounds(181, 11, 203, 37);
          panelProtocol.add(txtProtocolId);
        }
        {
          txtcProtocolId.setForeground(Color.DARK_GRAY);
          nimbusFix(Color.LIGHT_GRAY, txtcProtocolId);
          txtcProtocolId.setText("Protocol ID");
          txtcProtocolId.setFont(Application.TEXT_FONT);
          txtcProtocolId.setEditable(false);
          txtcProtocolId.setBounds(10, 11, 161, 37);
          panelProtocol.add(txtcProtocolId);
        }
      }
      {
        panelAlignment.setBounds(10, 81, 354, 92);
        panel1.add(panelAlignment);
        panelAlignment.setForeground(Color.DARK_GRAY);
        panelAlignment.setBackground(Color.LIGHT_GRAY);
        panelAlignment.setLayout(null);
        {
          txtcAlignmentAlgorithm.setForeground(Color.DARK_GRAY);
          nimbusFix(Color.LIGHT_GRAY, txtcAlignmentAlgorithm);
          txtcAlignmentAlgorithm.setText("Alignment Algorithm");
          txtcAlignmentAlgorithm.setEditable(false);
          txtcAlignmentAlgorithm.setBounds(10, 0, 334, 20);
          panelAlignment.add(txtcAlignmentAlgorithm);
        }
        comboAlignmentAlgo.setFont(Application.TEXT_FONT);
        comboAlignmentAlgo.setModel(new DefaultComboBoxModel<String>(new String[] { "BWA",
            "Bowtie", "Bowtie2", "Cushaw", "Cushaw-GPU", "TopHat" }));
        comboAlignmentAlgo.setBounds(10, 28, 334, 42);
        comboAlignmentAlgo.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent comboAction) {
            if (String.valueOf(comboAlignmentAlgo.getSelectedItem()).contains("Bowtie")) {
              checkRsem.setEnabled(true);
            } else {
              checkRsem.setSelected(false);
              checkRsem.setEnabled(false);
            }
          }
        });
        panelAlignment.add(comboAlignmentAlgo);
      }
      {
        panelCounting.setBounds(368, 81, 321, 92);
        panel1.add(panelCounting);
        panelCounting.setForeground(Color.DARK_GRAY);
        panelCounting.setBackground(Color.LIGHT_GRAY);
        panelCounting.setLayout(null);
        {
          txtcCountingMethods.setForeground(Color.DARK_GRAY);
          nimbusFix(Color.LIGHT_GRAY, txtcCountingMethods);
          txtcCountingMethods.setText("Counting Methods");
          txtcCountingMethods.setEditable(false);
          txtcCountingMethods.setBounds(10, 0, 301, 20);
          panelCounting.add(txtcCountingMethods);
        }
        checkFeatureCount.setFont(Application.TEXT_FONT);
        checkFeatureCount.setForeground(Color.DARK_GRAY);
        checkFeatureCount.setBackground(Color.LIGHT_GRAY);
        checkFeatureCount.setBounds(10, 27, 138, 23);

        panelCounting.add(checkFeatureCount);
        checkHtSeq.setFont(Application.TEXT_FONT);
        checkHtSeq.setForeground(Color.DARK_GRAY);
        checkHtSeq.setBackground(Color.LIGHT_GRAY);
        checkHtSeq.setBounds(10, 48, 138, 23);

        panelCounting.add(checkHtSeq);
        checkRsem.setFont(Application.TEXT_FONT);
        checkRsem.setForeground(Color.DARK_GRAY);
        checkRsem.setBackground(Color.LIGHT_GRAY);
        checkRsem.setBounds(10, 69, 138, 23);

        panelCounting.add(checkRsem);
      }
      {
        txtchPipelineOptions.setBounds(10, 11, 284, 42);
        panel1.add(txtchPipelineOptions);
        txtchPipelineOptions.setForeground(Color.DARK_GRAY);
        nimbusFix(Color.LIGHT_GRAY, txtchPipelineOptions);
        txtchPipelineOptions.setText("Pipeline Options");
        txtchPipelineOptions.setFont(Application.HEADER_FONT);
        txtchPipelineOptions.setEditable(false);
      }
      {
        panelNumberOfCores.setForeground(Color.DARK_GRAY);
        panelNumberOfCores.setBackground(Color.LIGHT_GRAY);
        panelNumberOfCores.setBounds(10, 277, 483, 40);
        panelManualHolder.add(panelNumberOfCores);
        panelNumberOfCores.setLayout(null);
        spinNumberOfCores.setModel(new SpinnerNumberModel(new Integer(4), new Integer(0), null,
            new Integer(1)));
        spinNumberOfCores.setBounds(380, 0, 77, 40);
        panelNumberOfCores.add(spinNumberOfCores);
        {
          txtcNumberOfCores.setForeground(Color.DARK_GRAY);
          nimbusFix(Color.LIGHT_GRAY, txtcNumberOfCores);
          txtcNumberOfCores.setEditable(false);
          txtcNumberOfCores.setBounds(0, 0, 370, 40);
          panelNumberOfCores.add(txtcNumberOfCores);
          txtcNumberOfCores.setText("Number of Cores to Use");
        }
      }
      {
        panelParallelExpressionComp.setForeground(Color.DARK_GRAY);
        panelParallelExpressionComp.setBackground(Color.LIGHT_GRAY);
        panelParallelExpressionComp.setLayout(null);
        panelParallelExpressionComp.setBounds(526, 277, 456, 40);
        panelManualHolder.add(panelParallelExpressionComp);
        {
          spinParallelExpression.setModel(new SpinnerNumberModel(new Integer(1), new Integer(0),
              null, new Integer(1)));
          spinParallelExpression.setBounds(380, 0, 76, 40);
          panelParallelExpressionComp.add(spinParallelExpression);
        }
        {
          txtcParallelExpressionComp.setForeground(Color.DARK_GRAY);
          nimbusFix(Color.LIGHT_GRAY, txtcParallelExpressionComp);
          txtcParallelExpressionComp.setEditable(false);
          txtcParallelExpressionComp
              .setText("Parallel Computation Streams for Expression Computation");
          txtcParallelExpressionComp.setBounds(0, 0, 370, 40);
          panelParallelExpressionComp.add(txtcParallelExpressionComp);
        }
      }
      {
        panelParallelDataComp.setForeground(Color.DARK_GRAY);
        panelParallelDataComp.setBackground(Color.LIGHT_GRAY);
        panelParallelDataComp.setLayout(null);
        panelParallelDataComp.setBounds(10, 328, 483, 40);
        panelManualHolder.add(panelParallelDataComp);
        {
          spinParallelDataPrep.setModel(new SpinnerNumberModel(new Integer(1), new Integer(0),
              null, new Integer(1)));
          spinParallelDataPrep.setBounds(380, 0, 77, 40);
          panelParallelDataComp.add(spinParallelDataPrep);
        }
        {
          txtParallelDataPrep.setForeground(Color.DARK_GRAY);
          nimbusFix(Color.LIGHT_GRAY, txtParallelDataPrep);
          txtParallelDataPrep.setEditable(false);
          txtParallelDataPrep.setText("Parallel Computation Streams for Data Preparation");
          txtParallelDataPrep.setBounds(0, 0, 370, 40);
          panelParallelDataComp.add(txtParallelDataPrep);
        }
      }
      {
        panelMaxFragLen.setForeground(Color.DARK_GRAY);
        panelMaxFragLen.setBackground(Color.LIGHT_GRAY);
        panelMaxFragLen.setLayout(null);
        panelMaxFragLen.setBounds(526, 328, 456, 40);
        panelManualHolder.add(panelMaxFragLen);
        {
          spinMaxFragLen.setModel(new SpinnerNumberModel(new Integer(1000), new Integer(0), null,
              new Integer(100)));
          spinMaxFragLen.setBounds(380, 0, 76, 40);
          panelMaxFragLen.add(spinMaxFragLen);
        }
        {
          txtcMaximalFragmentLength.setForeground(Color.DARK_GRAY);
          nimbusFix(Color.LIGHT_GRAY, txtcMaximalFragmentLength);
          txtcMaximalFragmentLength.setEditable(false);
          txtcMaximalFragmentLength.setText("Maximal Fragment Length");
          txtcMaximalFragmentLength.setBounds(0, 0, 370, 40);
          panelMaxFragLen.add(txtcMaximalFragmentLength);
        }
      }
      {
        panelMaxBufferConfInts.setForeground(Color.DARK_GRAY);
        panelMaxBufferConfInts.setBackground(Color.LIGHT_GRAY);
        panelMaxBufferConfInts.setLayout(null);
        panelMaxBufferConfInts.setBounds(10, 384, 483, 40);
        panelManualHolder.add(panelMaxBufferConfInts);
        {
          spinMaxBuffer.setModel(new SpinnerNumberModel(new Integer(4096), new Integer(0), null,
              new Integer(1024)));
          spinMaxBuffer.setBounds(379, 0, 78, 40);
          panelMaxBufferConfInts.add(spinMaxBuffer);
        }
        {
          txtcMaximumAuxiliaryBuffer.setForeground(Color.DARK_GRAY);
          nimbusFix(Color.LIGHT_GRAY, txtcMaximumAuxiliaryBuffer);
          txtcMaximumAuxiliaryBuffer.setEditable(false);
          txtcMaximumAuxiliaryBuffer
              .setText("Maximum Auxiliary Buffer for Computing Credibility Intervals");
          txtcMaximumAuxiliaryBuffer.setBounds(0, 0, 370, 40);
          panelMaxBufferConfInts.add(txtcMaximumAuxiliaryBuffer);
        }
      }
      {
        btnExtractCoverage.setFont(Application.TEXT_FONT);
        btnExtractCoverage.setBounds(10, 229, 983, 33);
        btnExtractCoverage.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent buttonAction) {
            if (btnExtractCoverage.getText().equals(ButtonEnums.OptionButton.EXTRACT.value)) {
              btnExtractCoverage.setText(ButtonEnums.OptionButton.NO_EXTRACT.value);
            } else {
              btnExtractCoverage.setText(ButtonEnums.OptionButton.EXTRACT.value);
            }
          }
        });
        panelManualHolder.add(btnExtractCoverage);
      }
      {
        btnOutputGenomeBam.setFont(Application.TEXT_FONT);
        btnOutputGenomeBam.setBounds(10, 178, 983, 40);
        btnOutputGenomeBam.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent buttonAction) {
            if (btnOutputGenomeBam.getText().equals(ButtonEnums.OptionButton.OUTPUT.value)) {
              btnOutputGenomeBam.setText(ButtonEnums.OptionButton.NO_OUTPUT.value);
            } else {
              btnOutputGenomeBam.setText(ButtonEnums.OptionButton.OUTPUT.value);
            }
          }
        });
        panelManualHolder.add(btnOutputGenomeBam);
      }
    }
    {
      buttonPane.setBackground(Color.GRAY);
      buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
      getContentPane().add(buttonPane, BorderLayout.SOUTH);
      {
        okButton.setActionCommand("OK");
        okButton.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent arg0) {
            Application.att.setaAlgor(String.valueOf(comboAlignmentAlgo.getSelectedItem()));
            Application.att.setRGenome(txtReferenceGenome.getText());
            Application.att.setFasta(txtReferenceFasta.getText());
            Application.att.setGffname(txtReferenceFeatures.getText());
            Application.att.setIdAttr(txtReferenceId.getText());
            Application.att.setFeatureColumn(String.valueOf(spinColumnGtf.getValue()));
            //
            if (String.valueOf(comboAlignmentAlgo.getSelectedItem()).equals("Cushaw-GPU")) {
              Application.att.setGpuAccel("TRUE");
            } else {
              Application.att.setGpuAccel("FALSE");
            }
            //
            Application.att.setCores(String.valueOf(spinNumberOfCores.getValue()));
            if (checkFeatureCount.isSelected()) {
              Application.att.setFeatureCounts("FeatureCounts");
            } else {
              Application.att.setFeatureCounts("");
            }
            //
            if (checkRsem.isSelected()) {
              Application.att.setRsem("RSEM");
            } else {
              Application.att.setRsem("");
            }
            //
            if (checkHtSeq.isSelected()) {
              Application.att.setHtseq("HTSeq");
            } else {
              Application.att.setHtseq("");
            }
            //
            Application.run.setProtocolId(txtProtocolId.getText());
            Application.att.setStreams(String.valueOf(spinParallelExpression.getValue()));
            Application.att.setStreamsDataPrep(String.valueOf(spinParallelDataPrep.getValue()));
            Application.att.setFragMaxLength(String.valueOf(spinMaxFragLen.getValue()));
            Application.att.setCiMem(String.valueOf(spinMaxBuffer.getValue()));
            //
            if (btnExtractCoverage.getText().equals(ButtonEnums.OptionButton.EXTRACT.value)) {
              Application.att.setStrandExtract("TRUE");
            } else {
              Application.att.setStrandExtract("FALSE");
            }
            //
            if (btnOutputGenomeBam.getText().equals(ButtonEnums.OptionButton.OUTPUT.value)) {
              Application.att.setGenoBam("TRUE");
            } else {
              Application.att.setGenoBam("FALSE");
            }
            //
            dispose();
          }
        });
        buttonPane.add(okButton);
        getRootPane().setDefaultButton(okButton);
      }
      {
        cancelButton.setActionCommand("Cancel");
        buttonPane.add(cancelButton);
        cancelButton.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent arg0) {
            dispose();
          }
        });
      }
    }
  }
  //
  // Fixes a bug in Nimbus where it overrides the desired JTextPane
  // Background colors
  //
  void nimbusFix(Color background, JTextPane pane) {
    UIDefaults defaults = new UIDefaults();
    defaults.put("TextPane[Enabled].backgroundPainter", background);
    pane.putClientProperty("Nimbus.Overrides", defaults);
    pane.putClientProperty("Nimbus.Overrides.InheritDefaults", true);
    pane.setBackground(background);
  }
}
