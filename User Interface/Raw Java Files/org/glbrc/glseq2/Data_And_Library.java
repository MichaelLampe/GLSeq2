package org.glbrc.glseq2;

import java.awt.FlowLayout;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.UIDefaults;
import javax.swing.JTextPane;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.ImageIcon;
import javax.swing.JTextArea;

import java.awt.Color;
import java.io.File;

import javax.swing.JComboBox;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

public class Data_And_Library extends JDialog {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private final JTextArea rawDirectory = new JTextArea(GLSeq2_Main_Application.att.getDirectory());
	private final JTextPane txtpnDirectoryContainingRaw = new JTextPane();
	private final JPanel panel_1 = new JPanel();
	private final JButton preprocFilesButton = new JButton("");
	private final JTextArea preprocessedFiles = new JTextArea(GLSeq2_Main_Application.att.getDirectoryFq());
	private final JTextPane txtpnDirectoryContainingPreprocessed = new JTextPane();
	private final JPanel panel_2 = new JPanel();
	private final JButton baseDirButton = new JButton("");
	private final JTextArea destinationDirectory = new JTextArea(GLSeq2_Main_Application.att.getDestinationDirectory());
	private final JTextPane txtpnBaseOfDestination = new JTextPane();
	private final JPanel panel_3 = new JPanel();
	private final JButton rawFilesButton = new JButton("");
	private final JTextArea rawFileNames = new JTextArea(GLSeq2_Main_Application.att.getRawFileNames());
	private final JTextPane txtpnRawFileNames = new JTextPane();
	private final JButton zippedButton = new JButton();
	private final JButton endButton = new JButton();
	private final JPanel panel_4 = new JPanel();
	private final JTextArea strain = new JTextArea(GLSeq2_Main_Application.att.getStrain());
	private final JTextPane txtpnStrain = new JTextPane();
	private final JPanel panel_5 = new JPanel();
	private final JTextArea subsetofLibraries = new JTextArea(GLSeq2_Main_Application.att.getLibList());
	private final JTextPane txtpnSubsetOfLibraries = new JTextPane();
	private final JPanel panel_6 = new JPanel();
	private final JPanel panel_7 = new JPanel();
	private final JPanel panel_8 = new JPanel();
	private final JPanel panel_9 = new JPanel();
	private final JTextPane txtpnQualityScoresFormat = new JTextPane();
	private final JTextPane txtpnSupportedSequencingPlatforms = new JTextPane();
	private final JTextPane txtpnLibraryIdLength = new JTextPane();
	private final JTextPane txtpnStrandednessOfThe = new JTextPane();
	private final JComboBox<String> qualityScores = new JComboBox<String>();
	private final JComboBox<String> sequencingPlatforms = new JComboBox<String>();
	private final JComboBox<String> strandedness = new JComboBox<String>();
	private final JSpinner libraryIdLength = new JSpinner();
	private final JPanel panel_10 = new JPanel();
	private final JButton countableSamDirButton = new JButton("");
	private final JTextArea countableSamDir = new JTextArea(GLSeq2_Main_Application.att.getCountableSamDir());
	private final JTextPane txtpnCountableSamFile = new JTextPane();

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		try {
			Data_And_Library dialog = new Data_And_Library();
			dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
			dialog.setVisible(true);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Create the dialog.
	 */
	public Data_And_Library() {
		initGUI();
		// Initialize some variables
		if (!GLSeq2_Main_Application.att.getQScores().equals("")) {
			qualityScores.setSelectedItem(GLSeq2_Main_Application.att.getQScores());
		}
		if (!GLSeq2_Main_Application.att.getSeqPlatform().equals("")) {
			sequencingPlatforms.setSelectedItem(GLSeq2_Main_Application.att.getSeqPlatform());
		}
		if (!GLSeq2_Main_Application.att.getStrandExtract().equals("")) {
			strandedness.setSelectedItem(GLSeq2_Main_Application.att.getLibstrand());
		}
		//
		//
		if (GLSeq2_Main_Application.att.getLibNchar().equals("0")){
			libraryIdLength.setValue(4);
		} else{
			libraryIdLength.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getLibNchar()));
		}
		//
		//
		if (GLSeq2_Main_Application.att.getUnzipped().equals("TRUE")){
			zippedButton.setText("Using Unzipped Files (.fq)");
		} else{
			zippedButton.setText("Using Zipped Files (.gz)");
		}
		if (GLSeq2_Main_Application.att.getPairedEnd().equals("FALSE")){
			endButton.setText("Using Single Ended Data");
		} else{
			endButton.setText("Using Paired Ended Data");
		}
	}
	private void initGUI() {
		getContentPane().setBackground(Color.LIGHT_GRAY);
		setResizable(false);
		setBounds(100, 100, 840, 629);
		getContentPane().setLayout(null);
		{
			JPanel buttonPane = new JPanel();
			buttonPane.setBackground(Color.GRAY);
			buttonPane.setBounds(0, 568, 834, 33);
			buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
			getContentPane().add(buttonPane);
			{
				final JButton okButton = new JButton("Apply and Close");
				okButton.setActionCommand("OK");
				buttonPane.add(okButton);
				okButton.addActionListener(new ActionListener(){
					public void actionPerformed(ActionEvent arg0) {
						GLSeq2_Main_Application.att.setDirectory(rawDirectory.getText());
						GLSeq2_Main_Application.att.setDirectoryFq(preprocessedFiles.getText());
						GLSeq2_Main_Application.att.setDestinationDirectory(destinationDirectory.getText());
						GLSeq2_Main_Application.att.setRawFileNames(rawFileNames.getText());	
						GLSeq2_Main_Application.att.setCountableSamDir(countableSamDir.getText());
						//
						//private final JButton zippedButton = new JButton("Using Zipped Files (.gz)");
						//private final JButton endButton = new JButton("Using Paired Ended Data");
						if (endButton.getText().equals("Using Paired Ended Data")){
							GLSeq2_Main_Application.att.setPairedEnd("TRUE");
						}else{
							GLSeq2_Main_Application.att.setPairedEnd("FALSE");
						}
						//
						if (zippedButton.getText().equals("Using Zipped Files (.gz)")){
							GLSeq2_Main_Application.att.setUnzipped("FALSE");
						} else{
							GLSeq2_Main_Application.att.setUnzipped("TRUE");
						}
						//GLSeq2_Main_Application.att.setUnzipped(unzipped);
						//GLSeq2_Main_Application.att.setPairedEnd(pairedEnd);
						//
						GLSeq2_Main_Application.att.setStrain(strain.getText());
						GLSeq2_Main_Application.att.setLibList(subsetofLibraries.getText());
						GLSeq2_Main_Application.att.setQScores(String.valueOf(qualityScores.getSelectedItem()));
						GLSeq2_Main_Application.att.setSeqPlatform(String.valueOf(sequencingPlatforms.getSelectedItem()));
						GLSeq2_Main_Application.att.setLibstrand(String.valueOf(strandedness.getSelectedItem()));
						GLSeq2_Main_Application.att.setLibNchar(String.valueOf(libraryIdLength.getValue()));					
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
			JPanel glowContainer = new JPanel();
			glowContainer.setBackground(Color.LIGHT_GRAY);
			glowContainer.setBounds(0, 0, 298, 569);
			getContentPane().add(glowContainer);
			glowContainer.setLayout(null);
		}
		{
			JPanel filledLibraryOptions = new JPanel();
			filledLibraryOptions.setBackground(Color.LIGHT_GRAY);
			filledLibraryOptions.setBounds(308, 0, 516, 569);
			getContentPane().add(filledLibraryOptions);
			filledLibraryOptions.setLayout(null);
			{
				JTextPane txtpnCurrentDataAnd = new JTextPane();
				txtpnCurrentDataAnd.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnCurrentDataAnd);
				txtpnCurrentDataAnd.setFont(new Font("Arial", Font.PLAIN, 20));
				txtpnCurrentDataAnd.setEditable(false);
				txtpnCurrentDataAnd.setText("Current Data and Library Options");
				txtpnCurrentDataAnd.setBounds(10, 5, 444, 31);
				filledLibraryOptions.add(txtpnCurrentDataAnd);
			}
			{
				JPanel panel = new JPanel();
				panel.setForeground(Color.DARK_GRAY);
				panel.setBackground(Color.LIGHT_GRAY);
				panel.setBounds(10, 42, 496, 42);
				filledLibraryOptions.add(panel);
				panel.setLayout(null);
				{
					JButton rawDirButton = new JButton("");
					rawDirButton.setIcon(new ImageIcon(Data_And_Library.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
					rawDirButton.setBounds(132, 0, 38, 37);
					rawDirButton.addActionListener(new ActionListener() {
						public void actionPerformed(ActionEvent arg0) {
							final JFileChooser chooser = new JFileChooser();
							chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
							chooser.showOpenDialog(rawDirectory);
							try {
								File file = chooser.getSelectedFile();
								rawDirectory.setText(file.getAbsolutePath());
							} catch (NullPointerException e) {
							}
						}
					});
					panel.add(rawDirButton);
				}
				rawDirectory.setBounds(180, 0, 306, 37);
				
				panel.add(rawDirectory);
				txtpnDirectoryContainingRaw.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnDirectoryContainingRaw);
				txtpnDirectoryContainingRaw.setEditable(false);
				txtpnDirectoryContainingRaw.setFont(new Font("Arial", Font.PLAIN, 11));
				txtpnDirectoryContainingRaw.setText("Directory Containing Raw Files");
				txtpnDirectoryContainingRaw.setBounds(0, 0, 112, 37);
				
				panel.add(txtpnDirectoryContainingRaw);
			}
			panel_1.setForeground(Color.DARK_GRAY);
			panel_1.setBackground(Color.LIGHT_GRAY);
			panel_1.setLayout(null);
			panel_1.setBounds(10, 91, 496, 42);
			
			filledLibraryOptions.add(panel_1);
			panel_2.setForeground(Color.DARK_GRAY);
			panel_2.setBackground(Color.LIGHT_GRAY);
			panel_2.setLayout(null);
			panel_2.setBounds(10, 144, 496, 42);
			
			filledLibraryOptions.add(panel_2);
			panel_3.setForeground(Color.DARK_GRAY);
			panel_3.setBackground(Color.LIGHT_GRAY);
			panel_3.setLayout(null);
			panel_3.setBounds(10, 197, 496, 42);
			
			filledLibraryOptions.add(panel_3);
			zippedButton.setFont(new Font("Arial", Font.PLAIN, 15));
			zippedButton.setBounds(10, 295, 226, 59);
			zippedButton.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent arg0) {
					if (zippedButton.getText().equals("Using Zipped Files (.gz)")){
						zippedButton.setText("Using Unzipped Files (.fq)");
					}
					else{
						zippedButton.setText("Using Zipped Files (.gz)");
					}
				}
			});
			
			filledLibraryOptions.add(zippedButton);
			endButton.setFont(new Font("Arial", Font.PLAIN, 15));
			endButton.setBounds(280, 295, 226, 59);
			endButton.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent arg0) {
					if (endButton.getText().equals("Using Paired Ended Data")){
						endButton.setText("Using Single Ended Data");
					}
					else{
						endButton.setText("Using Paired Ended Data");
					}
				}
			});
			
			filledLibraryOptions.add(endButton);
			panel_4.setForeground(Color.DARK_GRAY);
			panel_4.setBackground(Color.LIGHT_GRAY);
			panel_4.setLayout(null);
			panel_4.setBounds(10, 359, 496, 42);
			
			filledLibraryOptions.add(panel_4);
			panel_5.setForeground(Color.DARK_GRAY);
			panel_5.setBackground(Color.LIGHT_GRAY);
			panel_5.setLayout(null);
			panel_5.setBounds(10, 398, 496, 42);
			
			filledLibraryOptions.add(panel_5);
			panel_6.setForeground(Color.DARK_GRAY);
			panel_6.setBackground(Color.LIGHT_GRAY);
			panel_6.setBounds(10, 451, 240, 48);
			
			filledLibraryOptions.add(panel_6);
			panel_6.setLayout(null);
			txtpnQualityScoresFormat.setForeground(Color.DARK_GRAY);
			nimbusFix(Color.LIGHT_GRAY,txtpnQualityScoresFormat);
			txtpnQualityScoresFormat.setEditable(false);
			txtpnQualityScoresFormat.setFont(new Font("Arial", Font.PLAIN, 11));
			txtpnQualityScoresFormat.setText("Quality Scores Format");
			txtpnQualityScoresFormat.setBounds(10, 0, 220, 20);
			
			panel_6.add(txtpnQualityScoresFormat);
			qualityScores.setModel(new DefaultComboBoxModel<String>(new String[] {"phred33", "phred64"}));
			qualityScores.setBounds(10, 28, 220, 20);
			
			panel_6.add(qualityScores);
			panel_7.setForeground(Color.DARK_GRAY);
			panel_7.setBackground(Color.LIGHT_GRAY);
			panel_7.setBounds(266, 451, 240, 48);
			
			filledLibraryOptions.add(panel_7);
			panel_7.setLayout(null);
			txtpnSupportedSequencingPlatforms.setForeground(Color.DARK_GRAY);
			nimbusFix(Color.LIGHT_GRAY,txtpnSupportedSequencingPlatforms);
			txtpnSupportedSequencingPlatforms.setEditable(false);
			txtpnSupportedSequencingPlatforms.setText("Supported Sequencing Platforms");
			txtpnSupportedSequencingPlatforms.setFont(new Font("Arial", Font.PLAIN, 11));
			txtpnSupportedSequencingPlatforms.setBounds(10, 0, 220, 20);
			
			panel_7.add(txtpnSupportedSequencingPlatforms);
			sequencingPlatforms.setModel(new DefaultComboBoxModel<String>(new String[] {"illumina", "capillary", "ls454", "solid", "helicos", "iontorrent", "pacbio"}));
			sequencingPlatforms.setBounds(10, 28, 220, 20);
			
			panel_7.add(sequencingPlatforms);
			panel_8.setForeground(Color.DARK_GRAY);
			panel_8.setBackground(Color.LIGHT_GRAY);
			panel_8.setBounds(10, 510, 240, 48);
			
			filledLibraryOptions.add(panel_8);
			panel_8.setLayout(null);
			txtpnStrandednessOfThe.setForeground(Color.DARK_GRAY);
			nimbusFix(Color.LIGHT_GRAY,txtpnStrandednessOfThe);
			txtpnStrandednessOfThe.setEditable(false);
			txtpnStrandednessOfThe.setText("Strandedness of the Library");
			txtpnStrandednessOfThe.setFont(new Font("Arial", Font.PLAIN, 11));
			txtpnStrandednessOfThe.setBounds(10, 0, 220, 20);
			
			panel_8.add(txtpnStrandednessOfThe);
			strandedness.setModel(new DefaultComboBoxModel<String>(new String[] {"R", "F", "NULL"}));
			strandedness.setBounds(10, 28, 220, 20);
			
			panel_8.add(strandedness);
			panel_9.setForeground(Color.DARK_GRAY);
			panel_9.setBackground(Color.LIGHT_GRAY);
			panel_9.setBounds(266, 510, 240, 48);
			
			filledLibraryOptions.add(panel_9);
			panel_10.setLayout(null);
			panel_10.setForeground(Color.DARK_GRAY);
			panel_10.setBackground(Color.LIGHT_GRAY);
			panel_10.setBounds(10, 250, 496, 42);
			
			filledLibraryOptions.add(panel_10);
		}
		countableSamDirButton.setIcon(new ImageIcon(Data_And_Library.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
		countableSamDirButton.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				final JFileChooser chooser = new JFileChooser();
				chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				chooser.showOpenDialog(countableSamDir);
				try {
					File file = chooser.getSelectedFile();
					countableSamDir.setText(file.getAbsolutePath());
				} catch (NullPointerException ea) {
				}
			}
		});
		countableSamDirButton.setBounds(132, 0, 38, 37);
		
		panel_10.add(countableSamDirButton);
		countableSamDir.setBounds(180, 0, 306, 37);
		
		panel_10.add(countableSamDir);
		txtpnCountableSamFile.setText("Countable SAM File Directory");
		txtpnCountableSamFile.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtpnCountableSamFile);
		txtpnCountableSamFile.setFont(new Font("Arial", Font.PLAIN, 11));
		txtpnCountableSamFile.setEditable(false);
		txtpnCountableSamFile.setBounds(0, 0, 112, 37);
		
		panel_10.add(txtpnCountableSamFile);
		panel_9.setLayout(null);
		txtpnLibraryIdLength.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtpnLibraryIdLength);
		txtpnLibraryIdLength.setEditable(false);
		txtpnLibraryIdLength.setText("Library ID Length");
		txtpnLibraryIdLength.setFont(new Font("Arial", Font.PLAIN, 11));
		txtpnLibraryIdLength.setBounds(0, 11, 141, 37);
		
		panel_9.add(txtpnLibraryIdLength);
		libraryIdLength.setModel(new SpinnerNumberModel(new Integer(4), new Integer(0), null, new Integer(1)));
		libraryIdLength.setBounds(169, 0, 61, 48);
		
		panel_9.add(libraryIdLength);
		subsetofLibraries.setBounds(180, 11, 306, 31);
		
		panel_5.add(subsetofLibraries);
		txtpnSubsetOfLibraries.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtpnSubsetOfLibraries);
		txtpnSubsetOfLibraries.setEditable(false);
		txtpnSubsetOfLibraries.setText("Subset of Libraries to Process");
		txtpnSubsetOfLibraries.setFont(new Font("Arial", Font.PLAIN, 11));
		txtpnSubsetOfLibraries.setBounds(10, 11, 160, 31);
		
		panel_5.add(txtpnSubsetOfLibraries);
		strain.setBounds(180, 11, 306, 31);
		
		panel_4.add(strain);
		txtpnStrain.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtpnStrain);
		txtpnStrain.setEditable(false);
		txtpnStrain.setText("Strain");
		txtpnStrain.setFont(new Font("Arial", Font.PLAIN, 11));
		txtpnStrain.setBounds(10, 11, 160, 31);
		
		panel_4.add(txtpnStrain);
		rawFilesButton.setIcon(new ImageIcon(Data_And_Library.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
		rawFilesButton.setBounds(132, 0, 38, 37);
		rawFilesButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				final JFileChooser chooser = new JFileChooser();
				chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				chooser.showOpenDialog(rawFileNames);
				try {
					File file = chooser.getSelectedFile();
					rawFileNames.setText(file.getAbsolutePath());
				} catch (NullPointerException e) {
				}
			}
		});
		
		panel_3.add(rawFilesButton);
		rawFileNames.setBounds(180, 0, 306, 37);
		panel_3.add(rawFileNames);
		txtpnRawFileNames.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtpnRawFileNames);
		txtpnRawFileNames.setEditable(false);
		txtpnRawFileNames.setText("Raw File Names");
		txtpnRawFileNames.setFont(new Font("Arial", Font.PLAIN, 11));
		txtpnRawFileNames.setBounds(0, 0, 112, 37);
		
		panel_3.add(txtpnRawFileNames);
		baseDirButton.setIcon(new ImageIcon(Data_And_Library.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
		baseDirButton.setBounds(132, 0, 38, 37);
		baseDirButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				final JFileChooser chooser = new JFileChooser();
				chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				chooser.showOpenDialog(destinationDirectory);
				try {
					File file = chooser.getSelectedFile();
					destinationDirectory.setText(file.getAbsolutePath());
				} catch (NullPointerException e) {
				}
			}
		});
		
		panel_2.add(baseDirButton);
		destinationDirectory.setBounds(180, 0, 306, 37);
		
		panel_2.add(destinationDirectory);
		txtpnBaseOfDestination.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtpnBaseOfDestination);
		txtpnBaseOfDestination.setEditable(false);
		txtpnBaseOfDestination.setText("Base of Destination Directory");
		txtpnBaseOfDestination.setFont(new Font("Arial", Font.PLAIN, 11));
		txtpnBaseOfDestination.setBounds(0, 0, 112, 37);
		
		panel_2.add(txtpnBaseOfDestination);
		preprocFilesButton.setIcon(new ImageIcon(Data_And_Library.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
		preprocFilesButton.setBounds(132, 0, 38, 37);
		preprocFilesButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				final JFileChooser chooser = new JFileChooser();
				chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				chooser.showOpenDialog(preprocessedFiles);
				try {
					File file = chooser.getSelectedFile();
					preprocessedFiles.setText(file.getAbsolutePath());
				} catch (NullPointerException e) {
				}
			}
		});
		
		panel_1.add(preprocFilesButton);
		preprocessedFiles.setBounds(180, 0, 306, 37);
		
		panel_1.add(preprocessedFiles);
		txtpnDirectoryContainingPreprocessed.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtpnDirectoryContainingPreprocessed);
		txtpnDirectoryContainingPreprocessed.setEditable(false);
		txtpnDirectoryContainingPreprocessed.setText("Directory Containing Pre-Processed Files");
		txtpnDirectoryContainingPreprocessed.setFont(new Font("Arial", Font.PLAIN, 11));
		txtpnDirectoryContainingPreprocessed.setBounds(0, 0, 112, 37);
		
		panel_1.add(txtpnDirectoryContainingPreprocessed);
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
