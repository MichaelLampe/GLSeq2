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
	private final JSpinner spinner = new JSpinner();
	private final JTextArea textArea = new JTextArea();
	private final JComboBox comboBox = new JComboBox();
	private final JCheckBox chckbxNewCheckBox = new JCheckBox("Feature Counts");
	private final JCheckBox chckbxHtseq = new JCheckBox("HTSeq");
	private final JCheckBox chckbxRsem = new JCheckBox("RSEM");

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
					JButton button = new JButton("");
					button.setIcon(new ImageIcon(Script_Running_Options.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
					button.setBounds(132, 11, 38, 37);
					button.addActionListener(new ActionListener() {
						public void actionPerformed(ActionEvent arg0) {
							final JFileChooser chooser = new JFileChooser();
							chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
							chooser.showOpenDialog(textArea);
							try {
								File file = chooser.getSelectedFile();
								textArea.setText(file.getAbsolutePath());
							} catch (NullPointerException e) {
							}
						}
					});
					panel_1.add(button);
				}
				{
					textArea.setBounds(180, 11, 306, 37);
					panel_1.add(textArea);
				}
				{
					JTextPane txtpnGlseqScriptDirectory = new JTextPane();
					txtpnGlseqScriptDirectory.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY,txtpnGlseqScriptDirectory);
					txtpnGlseqScriptDirectory.setText("GLSeq Script Directory");
					txtpnGlseqScriptDirectory.setFont(new Font("Arial", Font.PLAIN, 11));
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
					JTextArea textArea = new JTextArea();
					textArea.setBounds(98, 11, 90, 37);
					panel_1.add(textArea);
				}
				{
					JTextPane txtpnProtocolId = new JTextPane();
					txtpnProtocolId.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY,txtpnProtocolId);
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
					nimbusFix(Color.LIGHT_GRAY,txtpnAlignmentAlgorithm);
					txtpnAlignmentAlgorithm.setText("Alignment Algorithm");
					txtpnAlignmentAlgorithm.setEditable(false);
					txtpnAlignmentAlgorithm.setBounds(10, 0, 220, 20);
					panel_1.add(txtpnAlignmentAlgorithm);
				}
				comboBox.setFont(new Font("Arial", Font.PLAIN, 11));
				comboBox.setModel(new DefaultComboBoxModel(new String[] {"BWA", "Bowtie", "Bowtie2", "Cushaw", "Cushaw-GPU"}));
				comboBox.setBounds(10, 28, 220, 42);
				
				panel_1.add(comboBox);
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
					nimbusFix(Color.LIGHT_GRAY,txtpnCountingMethods);
					txtpnCountingMethods.setText("Counting Methods");
					txtpnCountingMethods.setEditable(false);
					txtpnCountingMethods.setBounds(10, 0, 220, 20);
					panel_1.add(txtpnCountingMethods);
				}
				chckbxNewCheckBox.setFont(new Font("Arial", Font.PLAIN, 11));
				chckbxNewCheckBox.setForeground(Color.DARK_GRAY);
				chckbxNewCheckBox.setBackground(Color.LIGHT_GRAY);
				chckbxNewCheckBox.setBounds(10, 27, 122, 23);
				
				panel_1.add(chckbxNewCheckBox);
				chckbxHtseq.setFont(new Font("Arial", Font.PLAIN, 11));
				chckbxHtseq.setForeground(Color.DARK_GRAY);
				chckbxHtseq.setBackground(Color.LIGHT_GRAY);
				chckbxHtseq.setBounds(10, 48, 64, 23);
				
				panel_1.add(chckbxHtseq);
				chckbxRsem.setFont(new Font("Arial", Font.PLAIN, 11));
				chckbxRsem.setForeground(Color.DARK_GRAY);
				chckbxRsem.setBackground(Color.LIGHT_GRAY);
				chckbxRsem.setBounds(10, 69, 64, 23);
				
				panel_1.add(chckbxRsem);
			}
			{
				JTextPane txtpnPipelineOptions = new JTextPane();
				txtpnPipelineOptions.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnPipelineOptions);
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
				spinner.setModel(new SpinnerNumberModel(new Integer(4), new Integer(1), null, new Integer(1)));
				spinner.setBounds(310, 0, 98, 40);
				
				panel_1.add(spinner);
				{
					JTextPane txtpnNumberOfCores = new JTextPane();
					txtpnNumberOfCores.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY,txtpnNumberOfCores);
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
					JSpinner spinner_1 = new JSpinner();
					spinner_1.setModel(new SpinnerNumberModel(new Integer(1), new Integer(1), null, new Integer(1)));
					spinner_1.setBounds(310, 0, 98, 40);
					panel_1.add(spinner_1);
				}
				{
					JTextPane txtpnParallelComputationStreams = new JTextPane();
					txtpnParallelComputationStreams.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY,txtpnParallelComputationStreams);
					txtpnParallelComputationStreams.setEditable(false);
					txtpnParallelComputationStreams.setText("Parallel Computation Streams for Expression Computation");
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
					JSpinner spinner_1 = new JSpinner();
					spinner_1.setModel(new SpinnerNumberModel(new Integer(1), new Integer(1), null, new Integer(1)));
					spinner_1.setBounds(310, 0, 98, 40);
					panel_1.add(spinner_1);
				}
				{
					JTextPane txtpnParallelComputationStreams_1 = new JTextPane();
					txtpnParallelComputationStreams_1.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY,txtpnParallelComputationStreams_1);
					txtpnParallelComputationStreams_1.setEditable(false);
					txtpnParallelComputationStreams_1.setText("Parallel Computation Streams for Data Preparation");
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
					JSpinner spinner_1 = new JSpinner();
					spinner_1.setModel(new SpinnerNumberModel(new Integer(1000), new Integer(0), null, new Integer(100)));
					spinner_1.setBounds(310, 0, 98, 40);
					panel_1.add(spinner_1);
				}
				{
					JTextPane txtpnMaximalFragmentLength = new JTextPane();
					txtpnMaximalFragmentLength.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY,txtpnMaximalFragmentLength);
					txtpnMaximalFragmentLength.setEditable(false);
					txtpnMaximalFragmentLength.setText("Maximal Fragment Length");
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
					JSpinner spinner_1 = new JSpinner();
					spinner_1.setModel(new SpinnerNumberModel(new Integer(4096), new Integer(0), null, new Integer(1024)));
					spinner_1.setBounds(310, 0, 98, 40);
					panel_1.add(spinner_1);
				}
				{
					JTextPane txtpnMaximumAuxiliaryBuffer = new JTextPane();
					txtpnMaximumAuxiliaryBuffer.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY,txtpnMaximumAuxiliaryBuffer);
					txtpnMaximumAuxiliaryBuffer.setEditable(false);
					txtpnMaximumAuxiliaryBuffer.setText("Maximum Auxiliary Buffer for Computing Credibility Intervals");
					txtpnMaximumAuxiliaryBuffer.setBounds(0, 0, 307, 40);
					panel_1.add(txtpnMaximumAuxiliaryBuffer);
				}
			}
			{
				final JButton btnNewButton = new JButton("Extracting Forward and Reverse Coverage from Original BAM File");
				btnNewButton.setFont(new Font("Arial", Font.PLAIN, 11));
				btnNewButton.setBounds(10, 481, 496, 33);
				btnNewButton.addActionListener(new ActionListener(){
					public void actionPerformed(ActionEvent e) {
						if (btnNewButton.getText().equals("Extracting Forward and Reverse Coverage from Original BAM File")){
							btnNewButton.setText("Not Extracting Forward and Reverse Coverage from Original BAM File");
						} else{
							btnNewButton.setText("Extracting Forward and Reverse Coverage from Original BAM File");
						}
					}	
				});
				panel.add(btnNewButton);
			}
			{
				final JButton btnComputingConfidenceIntervals = new JButton("Computing Confidence Intervals");
				btnComputingConfidenceIntervals.setFont(new Font("Arial", Font.PLAIN, 11));
				btnComputingConfidenceIntervals.setBounds(10, 525, 240, 33);
				btnComputingConfidenceIntervals.addActionListener(new ActionListener(){
					public void actionPerformed(ActionEvent e) {
						if (btnComputingConfidenceIntervals.getText().equals("Computing Confidence Intervals")){
							btnComputingConfidenceIntervals.setText("Not Computing Confidence Intervals");
						} else{
							btnComputingConfidenceIntervals.setText("Computing Confidence Intervals");
						}
					}	
				});
				panel.add(btnComputingConfidenceIntervals);
			}
			{
				final JButton btnOutputtingGenomeBam = new JButton("Outputting Genome BAM");
				btnOutputtingGenomeBam.setFont(new Font("Arial", Font.PLAIN, 11));
				btnOutputtingGenomeBam.setBounds(266, 525, 240, 33);
				btnOutputtingGenomeBam.addActionListener(new ActionListener(){
					public void actionPerformed(ActionEvent e) {
						if (btnOutputtingGenomeBam.getText().equals("Outputting Genome BAM")){
							btnOutputtingGenomeBam.setText("Not Outputting Genome BAM");
						} else{
							btnOutputtingGenomeBam.setText("Outputting Genome BAM");
						}
					}	
				});
				panel.add(btnOutputtingGenomeBam);
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
	void nimbusFix(Color background,JTextPane pane){
		  UIDefaults defaults = new UIDefaults();
		  defaults.put("TextPane[Enabled].backgroundPainter", background);
		  pane.putClientProperty("Nimbus.Overrides", defaults);
		  pane.putClientProperty("Nimbus.Overrides.InheritDefaults", true);
		  pane.setBackground(background);
	}
}
