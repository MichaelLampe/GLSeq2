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

	private final JPanel contentPanel = new JPanel();
	//
	private final JSpinner numberOfCores = new JSpinner();
	private final JTextArea glseqDirectory = new JTextArea();
	private final JComboBox<String> alignmentAlgo = new JComboBox<String>();
	private final JCheckBox featureCountBox = new JCheckBox("Feature Counts");
	private final JCheckBox htseqBox = new JCheckBox("HTSeq");
	private final JCheckBox rsemBox = new JCheckBox("RSEM");
	private final JButton button = new JButton("");
	private final JTextArea protocolID = new JTextArea();
	private final JSpinner parallelStreamsExpression = new JSpinner();
	private final JSpinner parallelStreamsData = new JSpinner();
	private final JSpinner maxFragLength = new JSpinner();
	private final JSpinner maxBufferConfInt = new JSpinner();
	private final JButton extractCoverage = new JButton();
	private final JButton computeConfIntervals = new JButton();
	private final JButton outputGenomeBam = new JButton();
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
		if (!String.valueOf(alignmentAlgo.getSelectedItem()).contains("Bowtie")) {
			rsemBox.setSelected(false);
			rsemBox.setEnabled(false);
		}
		numberOfCores.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getCores()));
		glseqDirectory.setText(GLSeq2_Main_Application.att.getScriptDirectory());
		alignmentAlgo.setSelectedItem(GLSeq2_Main_Application.att.getaAlgor());
		if(GLSeq2_Main_Application.att.getFeatureCounts().contains("FeatureCounts")){
			featureCountBox.setSelected(true);
		}
		if(GLSeq2_Main_Application.att.getHTSeq().contains("HTSeq")){
			htseqBox.setSelected(true);
		}
		if(GLSeq2_Main_Application.att.getRSEM().contains("RSEM") && String.valueOf(alignmentAlgo.getSelectedItem()).contains("Bowtie")){
			rsemBox.setSelected(true);
		}
		protocolID.setText(GLSeq2_Main_Application.run.getProtocolId());
		parallelStreamsExpression.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getStreams()));
		parallelStreamsData.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getStreamsDataPrep()));
		maxFragLength.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getFragMaxLength()));
		maxBufferConfInt.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getCiMem()));
		
		if (GLSeq2_Main_Application.att.getStrandExtract().equals("TRUE")){
			extractCoverage.setText("Extracting Forward and Reverse Coverage from Original BAM File");
		} else{
			extractCoverage.setText("NOT Extracting Forward and Reverse Coverage from Original BAM File");
		}
		
		if (GLSeq2_Main_Application.att.getCompConf().equals("TRUE")){
			computeConfIntervals.setText("Computing Confidence Intervals");
		} else{
			computeConfIntervals.setText("NOT Computing Confidence Intervals");
		}
		
		if (GLSeq2_Main_Application.att.getGenoBam().equals("TRUE")){
			outputGenomeBam.setText("Outputting Genome BAM");
		} else{
			outputGenomeBam.setText("NOT Outputting Genome Bam");
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
			JPanel panel = new JPanel();
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setLayout(null);
			panel.setBounds(0, 0, 298, 569);
			contentPanel.add(panel);
		}
		{
			JPanel panel = new JPanel();
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setBounds(308, 0, 516, 569);
			contentPanel.add(panel);
			panel.setLayout(null);
			{
				JPanel panel_1 = new JPanel();
				panel_1.setForeground(Color.DARK_GRAY);
				panel_1.setBackground(Color.LIGHT_GRAY);
				panel_1.setLayout(null);
				panel_1.setBounds(10, 64, 496, 59);
				panel.add(panel_1);
				{
					button.setIcon(new ImageIcon(
							Script_Running_Options.class
									.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
					button.setBounds(132, 11, 38, 37);
					button.addActionListener(new ActionListener() {
						public void actionPerformed(ActionEvent arg0) {
							final JFileChooser chooser = new JFileChooser();
							chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
							chooser.showOpenDialog(glseqDirectory);
							try {
								File file = chooser.getSelectedFile();
								glseqDirectory.setText(file.getAbsolutePath());
							} catch (NullPointerException e) {
							}
						}
					});
					panel_1.add(button);
				}
				{
					glseqDirectory.setBounds(180, 11, 306, 37);
					panel_1.add(glseqDirectory);
				}
				{
					JTextPane txtpnGlseqScriptDirectory = new JTextPane();
					txtpnGlseqScriptDirectory.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtpnGlseqScriptDirectory);
					txtpnGlseqScriptDirectory.setText("GLSeq Script Directory");
					txtpnGlseqScriptDirectory.setFont(new Font("Arial",
							Font.PLAIN, 11));
					txtpnGlseqScriptDirectory.setEditable(false);
					txtpnGlseqScriptDirectory.setBounds(10, 11, 112, 37);
					panel_1.add(txtpnGlseqScriptDirectory);
				}
			}
			{
				JPanel panel_1 = new JPanel();
				panel_1.setForeground(Color.DARK_GRAY);
				panel_1.setBackground(Color.LIGHT_GRAY);
				panel_1.setLayout(null);
				panel_1.setBounds(301, 0, 205, 59);
				panel.add(panel_1);
				{
					protocolID.setBounds(98, 11, 90, 37);
					panel_1.add(protocolID);
				}
				{
					JTextPane txtpnProtocolId = new JTextPane();
					txtpnProtocolId.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtpnProtocolId);
					txtpnProtocolId.setText("Protocol ID");
					txtpnProtocolId.setFont(new Font("Arial", Font.PLAIN, 11));
					txtpnProtocolId.setEditable(false);
					txtpnProtocolId.setBounds(10, 11, 73, 37);
					panel_1.add(txtpnProtocolId);
				}
			}
			{
				JPanel panel_1 = new JPanel();
				panel_1.setForeground(Color.DARK_GRAY);
				panel_1.setBackground(Color.LIGHT_GRAY);
				panel_1.setBounds(10, 134, 240, 81);
				panel.add(panel_1);
				panel_1.setLayout(null);
				{
					JTextPane txtpnAlignmentAlgorithm = new JTextPane();
					txtpnAlignmentAlgorithm.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtpnAlignmentAlgorithm);
					txtpnAlignmentAlgorithm.setText("Alignment Algorithm");
					txtpnAlignmentAlgorithm.setEditable(false);
					txtpnAlignmentAlgorithm.setBounds(10, 0, 220, 20);
					panel_1.add(txtpnAlignmentAlgorithm);
				}
				alignmentAlgo.setFont(new Font("Arial", Font.PLAIN, 11));
				alignmentAlgo.setModel(new DefaultComboBoxModel<String>(
						new String[] { "BWA", "Bowtie", "Bowtie2", "Cushaw",
								"Cushaw-GPU" }));
				alignmentAlgo.setBounds(10, 28, 220, 42);
				alignmentAlgo.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						if (String.valueOf(alignmentAlgo.getSelectedItem())
								.contains("Bowtie")) {
							rsemBox.setEnabled(true);
						} else {
							rsemBox.setSelected(false);
							rsemBox.setEnabled(false);
						}
					}
				});
				panel_1.add(alignmentAlgo);
			}
			{
				JPanel panel_1 = new JPanel();
				panel_1.setForeground(Color.DARK_GRAY);
				panel_1.setBackground(Color.LIGHT_GRAY);
				panel_1.setLayout(null);
				panel_1.setBounds(266, 123, 240, 92);
				panel.add(panel_1);
				{
					JTextPane txtpnCountingMethods = new JTextPane();
					txtpnCountingMethods.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtpnCountingMethods);
					txtpnCountingMethods.setText("Counting Methods");
					txtpnCountingMethods.setEditable(false);
					txtpnCountingMethods.setBounds(10, 0, 220, 20);
					panel_1.add(txtpnCountingMethods);
				}
				featureCountBox.setFont(new Font("Arial", Font.PLAIN, 11));
				featureCountBox.setForeground(Color.DARK_GRAY);
				featureCountBox.setBackground(Color.LIGHT_GRAY);
				featureCountBox.setBounds(10, 27, 122, 23);

				panel_1.add(featureCountBox);
				htseqBox.setFont(new Font("Arial", Font.PLAIN, 11));
				htseqBox.setForeground(Color.DARK_GRAY);
				htseqBox.setBackground(Color.LIGHT_GRAY);
				htseqBox.setBounds(10, 48, 64, 23);

				panel_1.add(htseqBox);
				rsemBox.setFont(new Font("Arial", Font.PLAIN, 11));
				rsemBox.setForeground(Color.DARK_GRAY);
				rsemBox.setBackground(Color.LIGHT_GRAY);
				rsemBox.setBounds(10, 69, 64, 23);

				panel_1.add(rsemBox);
			}
			{
				JTextPane txtpnPipelineOptions = new JTextPane();
				txtpnPipelineOptions.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY, txtpnPipelineOptions);
				txtpnPipelineOptions.setText("Pipeline Options");
				txtpnPipelineOptions.setFont(new Font("Arial", Font.PLAIN, 20));
				txtpnPipelineOptions.setEditable(false);
				txtpnPipelineOptions.setBounds(10, 11, 284, 42);
				panel.add(txtpnPipelineOptions);
			}
			{
				JPanel panel_1 = new JPanel();
				panel_1.setForeground(Color.DARK_GRAY);
				panel_1.setBackground(Color.LIGHT_GRAY);
				panel_1.setBounds(10, 226, 496, 40);
				panel.add(panel_1);
				panel_1.setLayout(null);
				numberOfCores.setModel(new SpinnerNumberModel(new Integer(4), new Integer(0), null, new Integer(1)));
				numberOfCores.setBounds(310, 0, 98, 40);

				panel_1.add(numberOfCores);
				{
					JTextPane txtpnNumberOfCores = new JTextPane();
					txtpnNumberOfCores.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtpnNumberOfCores);
					txtpnNumberOfCores.setEditable(false);
					txtpnNumberOfCores.setBounds(0, 0, 300, 40);
					panel_1.add(txtpnNumberOfCores);
					txtpnNumberOfCores.setText("Number of Cores to Use");
				}
			}
			{
				JPanel panel_1 = new JPanel();
				panel_1.setForeground(Color.DARK_GRAY);
				panel_1.setBackground(Color.LIGHT_GRAY);
				panel_1.setLayout(null);
				panel_1.setBounds(10, 277, 496, 40);
				panel.add(panel_1);
				{
					parallelStreamsExpression.setModel(new SpinnerNumberModel(new Integer(1), new Integer(0), null, new Integer(1)));
					parallelStreamsExpression.setBounds(310, 0, 98, 40);
					panel_1.add(parallelStreamsExpression);
				}
				{
					JTextPane txtpnParallelComputationStreams = new JTextPane();
					txtpnParallelComputationStreams
							.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtpnParallelComputationStreams);
					txtpnParallelComputationStreams.setEditable(false);
					txtpnParallelComputationStreams
							.setText("Parallel Computation Streams for Expression Computation");
					txtpnParallelComputationStreams.setBounds(0, 0, 303, 40);
					panel_1.add(txtpnParallelComputationStreams);
				}
			}
			{
				JPanel panel_1 = new JPanel();
				panel_1.setForeground(Color.DARK_GRAY);
				panel_1.setBackground(Color.LIGHT_GRAY);
				panel_1.setLayout(null);
				panel_1.setBounds(10, 328, 496, 40);
				panel.add(panel_1);
				{
					parallelStreamsData.setModel(new SpinnerNumberModel(new Integer(1), new Integer(0), null, new Integer(1)));
					parallelStreamsData.setBounds(310, 0, 98, 40);
					panel_1.add(parallelStreamsData);
				}
				{
					JTextPane txtpnParallelComputationStreams_1 = new JTextPane();
					txtpnParallelComputationStreams_1
							.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY,
							txtpnParallelComputationStreams_1);
					txtpnParallelComputationStreams_1.setEditable(false);
					txtpnParallelComputationStreams_1
							.setText("Parallel Computation Streams for Data Preparation");
					txtpnParallelComputationStreams_1.setBounds(0, 0, 302, 40);
					panel_1.add(txtpnParallelComputationStreams_1);
				}
			}
			{
				JPanel panel_1 = new JPanel();
				panel_1.setForeground(Color.DARK_GRAY);
				panel_1.setBackground(Color.LIGHT_GRAY);
				panel_1.setLayout(null);
				panel_1.setBounds(10, 379, 496, 40);
				panel.add(panel_1);
				{
					maxFragLength.setModel(new SpinnerNumberModel(new Integer(1000), new Integer(0), null, new Integer(100)));
					maxFragLength.setBounds(310, 0, 98, 40);
					panel_1.add(maxFragLength);
				}
				{
					JTextPane txtpnMaximalFragmentLength = new JTextPane();
					txtpnMaximalFragmentLength.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtpnMaximalFragmentLength);
					txtpnMaximalFragmentLength.setEditable(false);
					txtpnMaximalFragmentLength
							.setText("Maximal Fragment Length");
					txtpnMaximalFragmentLength.setBounds(0, 0, 292, 40);
					panel_1.add(txtpnMaximalFragmentLength);
				}
			}
			{
				JPanel panel_1 = new JPanel();
				panel_1.setForeground(Color.DARK_GRAY);
				panel_1.setBackground(Color.LIGHT_GRAY);
				panel_1.setLayout(null);
				panel_1.setBounds(10, 430, 496, 40);
				panel.add(panel_1);
				{
					maxBufferConfInt.setModel(new SpinnerNumberModel(new Integer(4096), new Integer(0), null, new Integer(1024)));
					maxBufferConfInt.setBounds(310, 0, 98, 40);
					panel_1.add(maxBufferConfInt);
				}
				{
					JTextPane txtpnMaximumAuxiliaryBuffer = new JTextPane();
					txtpnMaximumAuxiliaryBuffer.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtpnMaximumAuxiliaryBuffer);
					txtpnMaximumAuxiliaryBuffer.setEditable(false);
					txtpnMaximumAuxiliaryBuffer
							.setText("Maximum Auxiliary Buffer for Computing Credibility Intervals");
					txtpnMaximumAuxiliaryBuffer.setBounds(0, 0, 307, 40);
					panel_1.add(txtpnMaximumAuxiliaryBuffer);
				}
			}
			{
				extractCoverage.setFont(new Font("Arial", Font.PLAIN, 11));
				extractCoverage.setBounds(10, 481, 496, 33);
				extractCoverage.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						if (extractCoverage
								.getText()
								.equals("Extracting Forward and Reverse Coverage from Original BAM File")) {
							extractCoverage
									.setText("Not Extracting Forward and Reverse Coverage from Original BAM File");
						} else {
							extractCoverage
									.setText("Extracting Forward and Reverse Coverage from Original BAM File");
						}
					}
				});
				panel.add(extractCoverage);
			}
			{
				computeConfIntervals.setFont(new Font("Arial", Font.PLAIN, 11));
				computeConfIntervals.setBounds(10, 525, 240, 33);
				computeConfIntervals.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						if (computeConfIntervals.getText().equals(
								"Computing Confidence Intervals")) {
							computeConfIntervals
									.setText("Not Computing Confidence Intervals");
						} else {
							computeConfIntervals
									.setText("Computing Confidence Intervals");
						}
					}
				});
				panel.add(computeConfIntervals);
			}
			{
				outputGenomeBam.setFont(new Font("Arial", Font.PLAIN, 11));
				outputGenomeBam.setBounds(266, 525, 240, 33);
				outputGenomeBam.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						if (outputGenomeBam.getText().equals(
								"Outputting Genome BAM")) {
							outputGenomeBam
									.setText("Not Outputting Genome BAM");
						} else {
							outputGenomeBam.setText("Outputting Genome BAM");
						}
					}
				});
				panel.add(outputGenomeBam);
			}
		}
		{
			JPanel buttonPane = new JPanel();
			buttonPane.setBackground(Color.GRAY);
			buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			{
				JButton okButton = new JButton("Apply and Close");
				okButton.setActionCommand("OK");
				okButton.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent arg0) {
						GLSeq2_Main_Application.att.setaAlgor(String.valueOf(alignmentAlgo.getSelectedItem()));
						GLSeq2_Main_Application.att.setCores(String.valueOf(numberOfCores.getValue()));
						GLSeq2_Main_Application.att.setScriptDirectory(glseqDirectory.getText());
						if (featureCountBox.isSelected()){
							GLSeq2_Main_Application.att.setFeatureCounts("FeatureCounts");
						} else{
							GLSeq2_Main_Application.att.setFeatureCounts("");
						}
						if (rsemBox.isSelected()){
							GLSeq2_Main_Application.att.setRSEM("RSEM");
						} else{
							GLSeq2_Main_Application.att.setRSEM("");
						}
						if (htseqBox.isSelected()){
							GLSeq2_Main_Application.att.setHTSeq("HTSeq");
						} else{
							GLSeq2_Main_Application.att.setHTSeq("");
						}
						GLSeq2_Main_Application.run.setProtocolId(protocolID.getText());
						GLSeq2_Main_Application.att.setStreams(String.valueOf(parallelStreamsExpression.getValue()));
						GLSeq2_Main_Application.att.setStreamsDataPrep(String.valueOf(parallelStreamsData.getValue()));
						GLSeq2_Main_Application.att.setFragMaxLength(String.valueOf(maxFragLength.getValue()));
						GLSeq2_Main_Application.att.setCiMem(String.valueOf(maxBufferConfInt.getValue()));
						
						if (extractCoverage.getText().equals("Extracting Forward and Reverse Coverage from Original BAM File")){
							GLSeq2_Main_Application.att.setStrandExtract("TRUE");
						} else{
							GLSeq2_Main_Application.att.setStrandExtract("FALSE");
						}
						
						if (computeConfIntervals.getText().equals("Computing Confidence Intervals")){
							GLSeq2_Main_Application.att.setCompConf("TRUE");
						} else{
							GLSeq2_Main_Application.att.setCompConf("FALSE");
						}
						
						if (outputGenomeBam.getText().equals("Outputting Genome BAM")){
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
