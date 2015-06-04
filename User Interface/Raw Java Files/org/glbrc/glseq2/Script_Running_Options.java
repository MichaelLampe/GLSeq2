package org.glbrc.glseq2;

import java.awt.BorderLayout;
import java.awt.FlowLayout;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.UIDefaults;
import javax.swing.border.EmptyBorder;
import javax.swing.JTextArea;
import javax.swing.JTextPane;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.ImageIcon;
import javax.swing.JSpinner;

import java.awt.Color;
import java.io.File;

import javax.swing.SpinnerNumberModel;
import javax.swing.JComboBox;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JCheckBox;

public class Script_Running_Options extends JDialog {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private final JPanel contentPanel = new JPanel();
	//
	private final JSpinner spinNumberOfCores = new JSpinner();
	private final JTextArea txtGlseqDirectory = new JTextArea();
	private final JComboBox<String> comboAlignmentAlgo = new JComboBox<String>();
	private final JCheckBox checkFeatureCount = new JCheckBox("Feature Counts");
	private final JCheckBox checkHtSeq = new JCheckBox("HTSeq");
	private final JCheckBox checkRsem = new JCheckBox("RSEM");
	private final JButton btnScriptDirectory = new JButton("");
	private final JTextArea txtProtocolID = new JTextArea();
	private final JSpinner spinParallelExpression = new JSpinner();
	private final JSpinner spinParallelDataPrep = new JSpinner();
	private final JSpinner spinMaxFragLen = new JSpinner();
	private final JSpinner spinMaxBuffer = new JSpinner();
	private final JButton btnExtractCoverage = new JButton();
	private final JButton btnComputeConfIntervals = new JButton();
	private final JButton btnOutputGenomeBam = new JButton();
	private final JPanel panelGlow = new JPanel();
	private final JPanel panelManualHolder = new JPanel();
	private final JPanel panelScript = new JPanel();
	private final JTextPane txtcGlseqScriptDirectory = new JTextPane();
	private final JPanel panelProtocol = new JPanel();
	private final JTextPane txtcProtocolId = new JTextPane();
	private final JPanel panelAlignment = new JPanel();
	private final JTextPane txtcAlignmentAlgorithm = new JTextPane();
	private final JPanel panelCounting = new JPanel();
	private final JTextPane txtcCountingMethods = new JTextPane();
	private final JTextPane txtchPipelineOptions = new JTextPane();
	private final JPanel panelNumberOfCores = new JPanel();
	private final JTextPane txtcNumberOfCores = new JTextPane();
	private final JPanel panelParallelExpressionComp = new JPanel();
	private final JTextPane txtcParallelComputationStreams = new JTextPane();
	private final JPanel panelParallelDataComp = new JPanel();
	private final JTextPane txtcParallelComputationStreams_1 = new JTextPane();
	private final JPanel panelMaxFragLen = new JPanel();
	private final JTextPane txtcMaximalFragmentLength = new JTextPane();
	private final JPanel panelMaxBufferConfInts = new JPanel();
	private final JTextPane txtcMaximumAuxiliaryBuffer = new JTextPane();
	private final JPanel buttonPane = new JPanel();
	private final JButton okButton = new JButton("Apply and Close");
	private final JButton cancelButton = new JButton("Cancel");
	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		try {
			Script_Running_Options dialog = new Script_Running_Options();
			dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
			dialog.setVisible(true);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Create the dialog.
	 */
	public Script_Running_Options() {
		initGUI();
		if (!String.valueOf(comboAlignmentAlgo.getSelectedItem()).contains("Bowtie")) {
			checkRsem.setSelected(false);
			checkRsem.setEnabled(false);
		}
		spinNumberOfCores.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getCores()));
		txtGlseqDirectory.setText(GLSeq2_Main_Application.att.getScriptDirectory());
		comboAlignmentAlgo.setSelectedItem(GLSeq2_Main_Application.att.getaAlgor());
		if(GLSeq2_Main_Application.att.getFeatureCounts().contains("FeatureCounts")){
			checkFeatureCount.setSelected(true);
		}
		if(GLSeq2_Main_Application.att.getHTSeq().contains("HTSeq")){
			checkHtSeq.setSelected(true);
		}
		if(GLSeq2_Main_Application.att.getRSEM().contains("RSEM") && String.valueOf(comboAlignmentAlgo.getSelectedItem()).contains("Bowtie")){
			checkRsem.setSelected(true);
		}
		txtProtocolID.setText(GLSeq2_Main_Application.run.getProtocolId());
		spinParallelExpression.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getStreams()));
		spinParallelDataPrep.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getStreamsDataPrep()));
		spinMaxFragLen.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getFragMaxLength()));
		spinMaxBuffer.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getCiMem()));
		
		if (GLSeq2_Main_Application.att.getStrandExtract().equals("TRUE")){
			btnExtractCoverage.setText(ButtonEnums.OptionButton.EXTRACT.value);
		} else{
			btnExtractCoverage.setText(ButtonEnums.OptionButton.NO_EXTRACT.value);
		}
		
		if (GLSeq2_Main_Application.att.getCompConf().equals("TRUE")){
			btnComputeConfIntervals.setText(ButtonEnums.OptionButton.COMPUTE.value);
		} else{
			btnComputeConfIntervals.setText(ButtonEnums.OptionButton.NO_COMPUTE.value);
		}
		
		if (GLSeq2_Main_Application.att.getGenoBam().equals("TRUE")){
			btnOutputGenomeBam.setText(ButtonEnums.OptionButton.OUTPUT.value);
		} else{
			btnOutputGenomeBam.setText(ButtonEnums.OptionButton.NO_OUTPUT.value);
		}
	}

	private void initGUI() {
		setResizable(false);
		setBounds(100, 100, 840, 629);
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBackground(Color.LIGHT_GRAY);
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		contentPanel.setLayout(null);
		{
			panelGlow.setBackground(Color.LIGHT_GRAY);
			panelGlow.setLayout(null);
			panelGlow.setBounds(0, 0, 298, 569);
			contentPanel.add(panelGlow);
		}
		{
			panelManualHolder.setBackground(Color.LIGHT_GRAY);
			panelManualHolder.setBounds(308, 0, 516, 569);
			contentPanel.add(panelManualHolder);
			panelManualHolder.setLayout(null);
			{
				panelScript.setForeground(Color.DARK_GRAY);
				panelScript.setBackground(Color.LIGHT_GRAY);
				panelScript.setLayout(null);
				panelScript.setBounds(10, 64, 496, 59);
				panelManualHolder.add(panelScript);
				{
					btnScriptDirectory.setIcon(new ImageIcon(
							Script_Running_Options.class
									.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
					btnScriptDirectory.setBounds(132, 11, 38, 37);
					btnScriptDirectory.addActionListener(new ActionListener() {
						public void actionPerformed(ActionEvent arg0) {
							final JFileChooser chooser = new JFileChooser();
							chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
							chooser.showOpenDialog(txtGlseqDirectory);
							try {
								File file = chooser.getSelectedFile();
								txtGlseqDirectory.setText(file.getAbsolutePath());
							} catch (NullPointerException e) {
							}
						}
					});
					panelScript.add(btnScriptDirectory);
				}
				{
					txtGlseqDirectory.setBounds(180, 11, 306, 37);
					panelScript.add(txtGlseqDirectory);
				}
				{
					txtcGlseqScriptDirectory.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtcGlseqScriptDirectory);
					txtcGlseqScriptDirectory.setText("GLSeq Script Directory");
					txtcGlseqScriptDirectory.setFont(new Font("Arial",
							Font.PLAIN, 11));
					txtcGlseqScriptDirectory.setEditable(false);
					txtcGlseqScriptDirectory.setBounds(10, 11, 112, 37);
					panelScript.add(txtcGlseqScriptDirectory);
				}
			}
			{
				panelProtocol.setForeground(Color.DARK_GRAY);
				panelProtocol.setBackground(Color.LIGHT_GRAY);
				panelProtocol.setLayout(null);
				panelProtocol.setBounds(301, 0, 205, 59);
				panelManualHolder.add(panelProtocol);
				{
					txtProtocolID.setBounds(98, 11, 90, 37);
					panelProtocol.add(txtProtocolID);
				}
				{
					txtcProtocolId.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtcProtocolId);
					txtcProtocolId.setText("Protocol ID");
					txtcProtocolId.setFont(GLSeq2_Main_Application.TEXT_FONT);
					txtcProtocolId.setEditable(false);
					txtcProtocolId.setBounds(10, 11, 73, 37);
					panelProtocol.add(txtcProtocolId);
				}
			}
			{
				panelAlignment.setForeground(Color.DARK_GRAY);
				panelAlignment.setBackground(Color.LIGHT_GRAY);
				panelAlignment.setBounds(10, 134, 240, 81);
				panelManualHolder.add(panelAlignment);
				panelAlignment.setLayout(null);
				{
					txtcAlignmentAlgorithm.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtcAlignmentAlgorithm);
					txtcAlignmentAlgorithm.setText("Alignment Algorithm");
					txtcAlignmentAlgorithm.setEditable(false);
					txtcAlignmentAlgorithm.setBounds(10, 0, 220, 20);
					panelAlignment.add(txtcAlignmentAlgorithm);
				}
				comboAlignmentAlgo.setFont(GLSeq2_Main_Application.TEXT_FONT);
				comboAlignmentAlgo.setModel(new DefaultComboBoxModel<String>(
						new String[] { "BWA", "Bowtie", "Bowtie2", "Cushaw",
								"Cushaw-GPU" }));
				comboAlignmentAlgo.setBounds(10, 28, 220, 42);
				comboAlignmentAlgo.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						if (String.valueOf(comboAlignmentAlgo.getSelectedItem())
								.contains("Bowtie")) {
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
				panelCounting.setForeground(Color.DARK_GRAY);
				panelCounting.setBackground(Color.LIGHT_GRAY);
				panelCounting.setLayout(null);
				panelCounting.setBounds(266, 123, 240, 92);
				panelManualHolder.add(panelCounting);
				{
					txtcCountingMethods.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtcCountingMethods);
					txtcCountingMethods.setText("Counting Methods");
					txtcCountingMethods.setEditable(false);
					txtcCountingMethods.setBounds(10, 0, 220, 20);
					panelCounting.add(txtcCountingMethods);
				}
				checkFeatureCount.setFont(GLSeq2_Main_Application.TEXT_FONT);
				checkFeatureCount.setForeground(Color.DARK_GRAY);
				checkFeatureCount.setBackground(Color.LIGHT_GRAY);
				checkFeatureCount.setBounds(10, 27, 122, 23);

				panelCounting.add(checkFeatureCount);
				checkHtSeq.setFont(GLSeq2_Main_Application.TEXT_FONT);
				checkHtSeq.setForeground(Color.DARK_GRAY);
				checkHtSeq.setBackground(Color.LIGHT_GRAY);
				checkHtSeq.setBounds(10, 48, 64, 23);

				panelCounting.add(checkHtSeq);
				checkRsem.setFont(GLSeq2_Main_Application.TEXT_FONT);
				checkRsem.setForeground(Color.DARK_GRAY);
				checkRsem.setBackground(Color.LIGHT_GRAY);
				checkRsem.setBounds(10, 69, 64, 23);

				panelCounting.add(checkRsem);
			}
			{
				txtchPipelineOptions.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY, txtchPipelineOptions);
				txtchPipelineOptions.setText("Pipeline Options");
				txtchPipelineOptions.setFont(GLSeq2_Main_Application.HEADER_FONT);
				txtchPipelineOptions.setEditable(false);
				txtchPipelineOptions.setBounds(10, 11, 284, 42);
				panelManualHolder.add(txtchPipelineOptions);
			}
			{
				panelNumberOfCores.setForeground(Color.DARK_GRAY);
				panelNumberOfCores.setBackground(Color.LIGHT_GRAY);
				panelNumberOfCores.setBounds(10, 226, 496, 40);
				panelManualHolder.add(panelNumberOfCores);
				panelNumberOfCores.setLayout(null);
				spinNumberOfCores.setModel(new SpinnerNumberModel(new Integer(4), new Integer(0), null, new Integer(1)));
				spinNumberOfCores.setBounds(310, 0, 98, 40);

				panelNumberOfCores.add(spinNumberOfCores);
				{
					txtcNumberOfCores.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtcNumberOfCores);
					txtcNumberOfCores.setEditable(false);
					txtcNumberOfCores.setBounds(0, 0, 300, 40);
					panelNumberOfCores.add(txtcNumberOfCores);
					txtcNumberOfCores.setText("Number of Cores to Use");
				}
			}
			{
				panelParallelExpressionComp.setForeground(Color.DARK_GRAY);
				panelParallelExpressionComp.setBackground(Color.LIGHT_GRAY);
				panelParallelExpressionComp.setLayout(null);
				panelParallelExpressionComp.setBounds(10, 277, 496, 40);
				panelManualHolder.add(panelParallelExpressionComp);
				{
					spinParallelExpression.setModel(new SpinnerNumberModel(new Integer(1), new Integer(0), null, new Integer(1)));
					spinParallelExpression.setBounds(310, 0, 98, 40);
					panelParallelExpressionComp.add(spinParallelExpression);
				}
				{
					txtcParallelComputationStreams
							.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtcParallelComputationStreams);
					txtcParallelComputationStreams.setEditable(false);
					txtcParallelComputationStreams
							.setText("Parallel Computation Streams for Expression Computation");
					txtcParallelComputationStreams.setBounds(0, 0, 303, 40);
					panelParallelExpressionComp.add(txtcParallelComputationStreams);
				}
			}
			{
				panelParallelDataComp.setForeground(Color.DARK_GRAY);
				panelParallelDataComp.setBackground(Color.LIGHT_GRAY);
				panelParallelDataComp.setLayout(null);
				panelParallelDataComp.setBounds(10, 328, 496, 40);
				panelManualHolder.add(panelParallelDataComp);
				{
					spinParallelDataPrep.setModel(new SpinnerNumberModel(new Integer(1), new Integer(0), null, new Integer(1)));
					spinParallelDataPrep.setBounds(310, 0, 98, 40);
					panelParallelDataComp.add(spinParallelDataPrep);
				}
				{
					txtcParallelComputationStreams_1
							.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY,
							txtcParallelComputationStreams_1);
					txtcParallelComputationStreams_1.setEditable(false);
					txtcParallelComputationStreams_1
							.setText("Parallel Computation Streams for Data Preparation");
					txtcParallelComputationStreams_1.setBounds(0, 0, 302, 40);
					panelParallelDataComp.add(txtcParallelComputationStreams_1);
				}
			}
			{
				panelMaxFragLen.setForeground(Color.DARK_GRAY);
				panelMaxFragLen.setBackground(Color.LIGHT_GRAY);
				panelMaxFragLen.setLayout(null);
				panelMaxFragLen.setBounds(10, 379, 496, 40);
				panelManualHolder.add(panelMaxFragLen);
				{
					spinMaxFragLen.setModel(new SpinnerNumberModel(new Integer(1000), new Integer(0), null, new Integer(100)));
					spinMaxFragLen.setBounds(310, 0, 98, 40);
					panelMaxFragLen.add(spinMaxFragLen);
				}
				{
					txtcMaximalFragmentLength.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtcMaximalFragmentLength);
					txtcMaximalFragmentLength.setEditable(false);
					txtcMaximalFragmentLength
							.setText("Maximal Fragment Length");
					txtcMaximalFragmentLength.setBounds(0, 0, 292, 40);
					panelMaxFragLen.add(txtcMaximalFragmentLength);
				}
			}
			{
				panelMaxBufferConfInts.setForeground(Color.DARK_GRAY);
				panelMaxBufferConfInts.setBackground(Color.LIGHT_GRAY);
				panelMaxBufferConfInts.setLayout(null);
				panelMaxBufferConfInts.setBounds(10, 430, 496, 40);
				panelManualHolder.add(panelMaxBufferConfInts);
				{
					spinMaxBuffer.setModel(new SpinnerNumberModel(new Integer(4096), new Integer(0), null, new Integer(1024)));
					spinMaxBuffer.setBounds(310, 0, 98, 40);
					panelMaxBufferConfInts.add(spinMaxBuffer);
				}
				{
					txtcMaximumAuxiliaryBuffer.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtcMaximumAuxiliaryBuffer);
					txtcMaximumAuxiliaryBuffer.setEditable(false);
					txtcMaximumAuxiliaryBuffer
							.setText("Maximum Auxiliary Buffer for Computing Credibility Intervals");
					txtcMaximumAuxiliaryBuffer.setBounds(0, 0, 307, 40);
					panelMaxBufferConfInts.add(txtcMaximumAuxiliaryBuffer);
				}
			}
			{
				btnExtractCoverage.setFont(GLSeq2_Main_Application.TEXT_FONT);
				btnExtractCoverage.setBounds(10, 481, 496, 33);
				btnExtractCoverage.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						if (btnExtractCoverage
								.getText()
								.equals(ButtonEnums.OptionButton.EXTRACT.value)) {
							btnExtractCoverage
									.setText(ButtonEnums.OptionButton.NO_EXTRACT.value);
						} else {
							btnExtractCoverage
									.setText(ButtonEnums.OptionButton.EXTRACT.value);
						}
					}
				});
				panelManualHolder.add(btnExtractCoverage);
			}
			{
				btnComputeConfIntervals.setFont(GLSeq2_Main_Application.TEXT_FONT);
				btnComputeConfIntervals.setBounds(10, 525, 240, 33);
				btnComputeConfIntervals.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						if (btnComputeConfIntervals.getText().equals(
								ButtonEnums.OptionButton.COMPUTE.value)) {
							btnComputeConfIntervals
									.setText(ButtonEnums.OptionButton.NO_COMPUTE.value);
						} else {
							btnComputeConfIntervals
									.setText(ButtonEnums.OptionButton.COMPUTE.value);
						}
					}
				});
				panelManualHolder.add(btnComputeConfIntervals);
			}
			{
				btnOutputGenomeBam.setFont(GLSeq2_Main_Application.TEXT_FONT);
				btnOutputGenomeBam.setBounds(266, 525, 240, 33);
				btnOutputGenomeBam.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						if (btnOutputGenomeBam.getText().equals(
								ButtonEnums.OptionButton.OUTPUT.value)) {
							btnOutputGenomeBam
									.setText(ButtonEnums.OptionButton.NO_OUTPUT.value);
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
						GLSeq2_Main_Application.att.setaAlgor(String.valueOf(comboAlignmentAlgo.getSelectedItem()));
						GLSeq2_Main_Application.att.setCores(String.valueOf(spinNumberOfCores.getValue()));
						GLSeq2_Main_Application.att.setScriptDirectory(txtGlseqDirectory.getText());
						if (checkFeatureCount.isSelected()){
							GLSeq2_Main_Application.att.setFeatureCounts("FeatureCounts");
						} else{
							GLSeq2_Main_Application.att.setFeatureCounts("");
						}
						if (checkRsem.isSelected()){
							GLSeq2_Main_Application.att.setRSEM("RSEM");
						} else{
							GLSeq2_Main_Application.att.setRSEM("");
						}
						if (checkHtSeq.isSelected()){
							GLSeq2_Main_Application.att.setHTSeq("HTSeq");
						} else{
							GLSeq2_Main_Application.att.setHTSeq("");
						}
						GLSeq2_Main_Application.run.setProtocolId(txtProtocolID.getText());
						GLSeq2_Main_Application.att.setStreams(String.valueOf(spinParallelExpression.getValue()));
						GLSeq2_Main_Application.att.setStreamsDataPrep(String.valueOf(spinParallelDataPrep.getValue()));
						GLSeq2_Main_Application.att.setFragMaxLength(String.valueOf(spinMaxFragLen.getValue()));
						GLSeq2_Main_Application.att.setCiMem(String.valueOf(spinMaxBuffer.getValue()));
						
						if (btnExtractCoverage.getText().equals(ButtonEnums.OptionButton.EXTRACT.value)){
							GLSeq2_Main_Application.att.setStrandExtract("TRUE");
						} else{
							GLSeq2_Main_Application.att.setStrandExtract("FALSE");
						}
						
						if (btnComputeConfIntervals.getText().equals(ButtonEnums.OptionButton.COMPUTE.value)){
							GLSeq2_Main_Application.att.setCompConf("TRUE");
						} else{
							GLSeq2_Main_Application.att.setCompConf("FALSE");
						}
						
						if (btnOutputGenomeBam.getText().equals(ButtonEnums.OptionButton.OUTPUT.value)){
							GLSeq2_Main_Application.att.setGenoBam("TRUE");
						} else{
							GLSeq2_Main_Application.att.setGenoBam("FALSE");
						}
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