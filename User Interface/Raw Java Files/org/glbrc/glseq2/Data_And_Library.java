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
	/**
	 * 
	 */
	private final JTextPane txtchCurrentDataAnd = new JTextPane();
	private final JTextArea txtrawDirectory = new JTextArea(GLSeq2_Main_Application.att.getDirectory());
	private final JTextPane txtcDirectoryContainingRaw = new JTextPane();
	private final JPanel panelPreProcessedFiles = new JPanel();
	private final JButton btnPreProcFiles = new JButton("");
	private final JTextArea txtPreProcessedFiles = new JTextArea(GLSeq2_Main_Application.att.getDirectoryFq());
	private final JTextPane txtcDirectoryContainingPreprocessed = new JTextPane();
	private final JPanel panelDestinationDir = new JPanel();
	private final JButton btnBaseDir = new JButton("");
	private final JTextArea txtDestinationDirectory = new JTextArea(GLSeq2_Main_Application.att.getDestinationDirectory());
	private final JTextPane txtcBaseOfDestination = new JTextPane();
	private final JPanel panelRawFiles = new JPanel();
	private final JButton btnRawFiles = new JButton("");
	private final JTextArea txtRawFileNames = new JTextArea(GLSeq2_Main_Application.att.getRawFileNames());
	private final JTextPane txtcRawFileNames = new JTextPane();
	private final JButton btnZipped = new JButton();
	private final JButton btnEnded = new JButton();
	private final JPanel panelStrain = new JPanel();
	private final JTextArea txtStrain = new JTextArea(GLSeq2_Main_Application.att.getStrain());
	private final JTextPane txtcStrain = new JTextPane();
	private final JPanel panelSubsetOfLibraries = new JPanel();
	private final JTextArea txtSubsetOfLibraries = new JTextArea(GLSeq2_Main_Application.att.getLibList());
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
	private final JTextArea txtCountableSamDirectory = new JTextArea(GLSeq2_Main_Application.att.getCountableSamDir());
	private final JTextPane txtcCountableSamFile = new JTextPane();

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
			comboQualityScores.setSelectedItem(GLSeq2_Main_Application.att.getQScores());
		}
		if (!GLSeq2_Main_Application.att.getSeqPlatform().equals("")) {
			comboSequencingPlatforms.setSelectedItem(GLSeq2_Main_Application.att.getSeqPlatform());
		}
		if (!GLSeq2_Main_Application.att.getStrandExtract().equals("")) {
			comboStrandedness.setSelectedItem(GLSeq2_Main_Application.att.getLibstrand());
		}
		//
		//
		if (GLSeq2_Main_Application.att.getLibNchar().equals("0")){
			spinLibraryIdLen.setValue(4);
		} else{
			spinLibraryIdLen.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getLibNchar()));
		}
		//
		//
		if (GLSeq2_Main_Application.att.getUnzipped().equals("TRUE")){
			btnZipped.setText(ButtonEnums.OptionButton.UNZIPPED.value);
		} else{
			btnZipped.setText(ButtonEnums.OptionButton.ZIPPED.value);
		}
		if (GLSeq2_Main_Application.att.getPairedEnd().equals("FALSE")){
			btnEnded.setText(ButtonEnums.OptionButton.SINGLE.value);
		} else{
			btnEnded.setText(ButtonEnums.OptionButton.PAIRED.value);
		}
	}
	private void initGUI() {
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
				okButton.addActionListener(new ActionListener(){
					public void actionPerformed(ActionEvent arg0) {
						/*
						 * Assigns all the values to the attribute file's fields
						 */
						GLSeq2_Main_Application.att.setDirectory(txtrawDirectory.getText());
						GLSeq2_Main_Application.att.setDirectoryFq(txtPreProcessedFiles.getText());
						GLSeq2_Main_Application.att.setDestinationDirectory(txtDestinationDirectory.getText());
						GLSeq2_Main_Application.att.setRawFileNames(txtRawFileNames.getText());	
						GLSeq2_Main_Application.att.setCountableSamDir(txtCountableSamDirectory.getText());
						//
						//private final JButton zippedButton = new JButton("Using Zipped Files (.gz)");
						//private final JButton endButton = new JButton("Using Paired Ended Data");
						if (btnEnded.getText().equals(ButtonEnums.OptionButton.PAIRED.value)){
							GLSeq2_Main_Application.att.setPairedEnd("TRUE");
						}else{
							GLSeq2_Main_Application.att.setPairedEnd("FALSE");
						}
						//
						if (btnZipped.getText().equals(ButtonEnums.OptionButton.ZIPPED.value)){
							GLSeq2_Main_Application.att.setUnzipped("FALSE");
						} else{
							GLSeq2_Main_Application.att.setUnzipped("TRUE");
						}
						//GLSeq2_Main_Application.att.setUnzipped(unzipped);
						//GLSeq2_Main_Application.att.setPairedEnd(pairedEnd);
						//
						GLSeq2_Main_Application.att.setStrain(txtStrain.getText());
						GLSeq2_Main_Application.att.setLibList(txtSubsetOfLibraries.getText());
						GLSeq2_Main_Application.att.setQScores(String.valueOf(comboQualityScores.getSelectedItem()));
						GLSeq2_Main_Application.att.setSeqPlatform(String.valueOf(comboSequencingPlatforms.getSelectedItem()));
						GLSeq2_Main_Application.att.setLibstrand(String.valueOf(comboStrandedness.getSelectedItem()));
						GLSeq2_Main_Application.att.setLibNchar(String.valueOf(spinLibraryIdLen.getValue()));					
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
			JPanel panelGlow = new JPanel();
			panelGlow.setBackground(Color.LIGHT_GRAY);
			panelGlow.setBounds(0, 0, 298, 569);
			getContentPane().add(panelGlow);
			panelGlow.setLayout(null);
		}
		{
			JPanel panelZipAndStrain = new JPanel();
			panelZipAndStrain.setBackground(Color.LIGHT_GRAY);
			panelZipAndStrain.setBounds(308, 0, 729, 587);
			getContentPane().add(panelZipAndStrain);
			panelZipAndStrain.setLayout(null);
			{
				txtchCurrentDataAnd.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtchCurrentDataAnd);
				txtchCurrentDataAnd.setFont(GLSeq2_Main_Application.HEADER_FONT);
				txtchCurrentDataAnd.setEditable(false);
				txtchCurrentDataAnd.setText("Current Data and Library Options");
				txtchCurrentDataAnd.setBounds(10, 5, 444, 48);
				panelZipAndStrain.add(txtchCurrentDataAnd);
			}
			{
				JPanel panelRawDirectory = new JPanel();
				panelRawDirectory.setForeground(Color.DARK_GRAY);
				panelRawDirectory.setBackground(Color.LIGHT_GRAY);
				panelRawDirectory.setBounds(10, 60, 719, 42);
				panelZipAndStrain.add(panelRawDirectory);
				panelRawDirectory.setLayout(null);
				{
					JButton btnRawDir = new JButton("");
					btnRawDir.setIcon(new ImageIcon(Data_And_Library.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
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
							}
						}
					});
					panelRawDirectory.add(btnRawDir);
				}
				txtrawDirectory.setBounds(319, 0, 400, 37);
				
				panelRawDirectory.add(txtrawDirectory);
				txtcDirectoryContainingRaw.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcDirectoryContainingRaw);
				txtcDirectoryContainingRaw.setEditable(false);
				txtcDirectoryContainingRaw.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcDirectoryContainingRaw.setText("Directory Containing Raw Files");
				txtcDirectoryContainingRaw.setBounds(0, 0, 261, 37);
				
				panelRawDirectory.add(txtcDirectoryContainingRaw);
			}
			panelPreProcessedFiles.setForeground(Color.DARK_GRAY);
			panelPreProcessedFiles.setBackground(Color.LIGHT_GRAY);
			panelPreProcessedFiles.setLayout(null);
			panelPreProcessedFiles.setBounds(10, 109, 719, 42);
			
			panelZipAndStrain.add(panelPreProcessedFiles);
			panelDestinationDir.setForeground(Color.DARK_GRAY);
			panelDestinationDir.setBackground(Color.LIGHT_GRAY);
			panelDestinationDir.setLayout(null);
			panelDestinationDir.setBounds(10, 162, 719, 42);
			
			panelZipAndStrain.add(panelDestinationDir);
			panelRawFiles.setForeground(Color.DARK_GRAY);
			panelRawFiles.setBackground(Color.LIGHT_GRAY);
			panelRawFiles.setLayout(null);
			panelRawFiles.setBounds(10, 215, 719, 42);
			
			panelZipAndStrain.add(panelRawFiles);
			btnZipped.setFont(new Font("Arial", Font.PLAIN, 15));
			btnZipped.setBounds(10, 313, 354, 59);
			btnZipped.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent arg0) {
					if (btnZipped.getText().equals(ButtonEnums.OptionButton.ZIPPED.value)){
						btnZipped.setText(ButtonEnums.OptionButton.UNZIPPED.value);
					}
					else{
						btnZipped.setText(ButtonEnums.OptionButton.ZIPPED.value);
					}
				}
			});
			
			panelZipAndStrain.add(btnZipped);
			btnEnded.setFont(new Font("Arial", Font.PLAIN, 15));
			btnEnded.setBounds(375, 313, 354, 59);
			btnEnded.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent arg0) {
					if (btnEnded.getText().equals(ButtonEnums.OptionButton.PAIRED.value)){
						btnEnded.setText(ButtonEnums.OptionButton.SINGLE.value);
					}
					else{
						btnEnded.setText(ButtonEnums.OptionButton.PAIRED.value);
					}
				}
			});
			
			panelZipAndStrain.add(btnEnded);
			panelStrain.setForeground(Color.DARK_GRAY);
			panelStrain.setBackground(Color.LIGHT_GRAY);
			panelStrain.setLayout(null);
			panelStrain.setBounds(10, 377, 719, 42);
			
			panelZipAndStrain.add(panelStrain);
			panelSubsetOfLibraries.setForeground(Color.DARK_GRAY);
			panelSubsetOfLibraries.setBackground(Color.LIGHT_GRAY);
			panelSubsetOfLibraries.setLayout(null);
			panelSubsetOfLibraries.setBounds(10, 416, 719, 42);
			
			panelZipAndStrain.add(panelSubsetOfLibraries);
			panelQualityScores.setForeground(Color.DARK_GRAY);
			panelQualityScores.setBackground(Color.LIGHT_GRAY);
			panelQualityScores.setBounds(10, 469, 354, 48);
			
			panelZipAndStrain.add(panelQualityScores);
			panelQualityScores.setLayout(null);
			txtcQualityScores.setForeground(Color.DARK_GRAY);
			nimbusFix(Color.LIGHT_GRAY,txtcQualityScores);
			txtcQualityScores.setEditable(false);
			txtcQualityScores.setFont(GLSeq2_Main_Application.TEXT_FONT);
			txtcQualityScores.setText("Quality Scores Format");
			txtcQualityScores.setBounds(10, 0, 334, 20);
			
			panelQualityScores.add(txtcQualityScores);
			comboQualityScores.setModel(new DefaultComboBoxModel<String>(new String[] {"phred33", "phred64"}));
			comboQualityScores.setBounds(10, 28, 334, 20);
			
			panelQualityScores.add(comboQualityScores);
			panelSequencing.setForeground(Color.DARK_GRAY);
			panelSequencing.setBackground(Color.LIGHT_GRAY);
			panelSequencing.setBounds(375, 469, 354, 48);
			
			panelZipAndStrain.add(panelSequencing);
			panelSequencing.setLayout(null);
			txtcSequencingPlatforms.setForeground(Color.DARK_GRAY);
			nimbusFix(Color.LIGHT_GRAY,txtcSequencingPlatforms);
			txtcSequencingPlatforms.setEditable(false);
			txtcSequencingPlatforms.setText("Supported Sequencing Platforms");
			txtcSequencingPlatforms.setFont(GLSeq2_Main_Application.TEXT_FONT);
			txtcSequencingPlatforms.setBounds(10, 0, 334, 20);
			
			panelSequencing.add(txtcSequencingPlatforms);
			comboSequencingPlatforms.setModel(new DefaultComboBoxModel<String>(new String[] {"illumina", "capillary", "ls454", "solid", "helicos", "iontorrent", "pacbio"}));
			comboSequencingPlatforms.setBounds(10, 28, 334, 20);
			
			panelSequencing.add(comboSequencingPlatforms);
			panelStrandedness.setForeground(Color.DARK_GRAY);
			panelStrandedness.setBackground(Color.LIGHT_GRAY);
			panelStrandedness.setBounds(10, 528, 354, 48);
			
			panelZipAndStrain.add(panelStrandedness);
			panelStrandedness.setLayout(null);
			txtcStrandedness.setForeground(Color.DARK_GRAY);
			nimbusFix(Color.LIGHT_GRAY,txtcStrandedness);
			txtcStrandedness.setEditable(false);
			txtcStrandedness.setText("Strandedness of the Library");
			txtcStrandedness.setFont(GLSeq2_Main_Application.TEXT_FONT);
			txtcStrandedness.setBounds(10, 0, 334, 20);
			
			panelStrandedness.add(txtcStrandedness);
			comboStrandedness.setModel(new DefaultComboBoxModel<String>(new String[] {"R", "F", "NULL"}));
			comboStrandedness.setBounds(10, 28, 334, 20);
			
			panelStrandedness.add(comboStrandedness);
			panelLibraryIdLength.setForeground(Color.DARK_GRAY);
			panelLibraryIdLength.setBackground(Color.LIGHT_GRAY);
			panelLibraryIdLength.setBounds(385, 528, 344, 48);
			
			panelZipAndStrain.add(panelLibraryIdLength);
			panelCountableSam.setLayout(null);
			panelCountableSam.setForeground(Color.DARK_GRAY);
			panelCountableSam.setBackground(Color.LIGHT_GRAY);
			panelCountableSam.setBounds(10, 268, 719, 42);
			
			panelZipAndStrain.add(panelCountableSam);
		}
		btnCountableSam.setIcon(new ImageIcon(Data_And_Library.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
		btnCountableSam.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				final JFileChooser chooser = new JFileChooser();
				chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				chooser.showOpenDialog(txtCountableSamDirectory);
				try {
					File file = chooser.getSelectedFile();
					txtCountableSamDirectory.setText(file.getAbsolutePath());
				} catch (NullPointerException ea) {
				}
			}
		});
		btnCountableSam.setBounds(271, 0, 38, 37);
		
		panelCountableSam.add(btnCountableSam);
		txtCountableSamDirectory.setBounds(319, 0, 400, 37);
		
		panelCountableSam.add(txtCountableSamDirectory);
		txtcCountableSamFile.setText("Countable SAM File Directory");
		txtcCountableSamFile.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtcCountableSamFile);
		txtcCountableSamFile.setFont(GLSeq2_Main_Application.TEXT_FONT);
		txtcCountableSamFile.setEditable(false);
		txtcCountableSamFile.setBounds(0, 0, 261, 37);
		
		panelCountableSam.add(txtcCountableSamFile);
		panelLibraryIdLength.setLayout(null);
		txtcLibraryIdLength.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtcLibraryIdLength);
		txtcLibraryIdLength.setEditable(false);
		txtcLibraryIdLength.setText("Library ID Length");
		txtcLibraryIdLength.setFont(GLSeq2_Main_Application.TEXT_FONT);
		txtcLibraryIdLength.setBounds(0, 11, 141, 37);
		
		panelLibraryIdLength.add(txtcLibraryIdLength);
		spinLibraryIdLen.setModel(new SpinnerNumberModel(new Integer(4), new Integer(0), null, new Integer(1)));
		spinLibraryIdLen.setBounds(169, 0, 61, 48);
		
		panelLibraryIdLength.add(spinLibraryIdLen);
		txtSubsetOfLibraries.setBounds(319, 11, 400, 31);
		
		panelSubsetOfLibraries.add(txtSubsetOfLibraries);
		txtcSubsetOfLibraries.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtcSubsetOfLibraries);
		txtcSubsetOfLibraries.setEditable(false);
		txtcSubsetOfLibraries.setText("Subset of Libraries to Process");
		txtcSubsetOfLibraries.setFont(GLSeq2_Main_Application.TEXT_FONT);
		txtcSubsetOfLibraries.setBounds(10, 11, 299, 31);
		
		panelSubsetOfLibraries.add(txtcSubsetOfLibraries);
		txtStrain.setBounds(319, 11, 400, 31);
		
		panelStrain.add(txtStrain);
		txtcStrain.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtcStrain);
		txtcStrain.setEditable(false);
		txtcStrain.setText("Strain");
		txtcStrain.setFont(GLSeq2_Main_Application.TEXT_FONT);
		txtcStrain.setBounds(10, 11, 299, 31);
		
		panelStrain.add(txtcStrain);
		btnRawFiles.setIcon(new ImageIcon(Data_And_Library.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
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
				}
			}
		});
		
		panelRawFiles.add(btnRawFiles);
		txtRawFileNames.setBounds(319, 0, 400, 37);
		panelRawFiles.add(txtRawFileNames);
		txtcRawFileNames.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtcRawFileNames);
		txtcRawFileNames.setEditable(false);
		txtcRawFileNames.setText("Raw File Names");
		txtcRawFileNames.setFont(GLSeq2_Main_Application.TEXT_FONT);
		txtcRawFileNames.setBounds(0, 0, 261, 37);
		
		panelRawFiles.add(txtcRawFileNames);
		btnBaseDir.setIcon(new ImageIcon(Data_And_Library.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
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
				}
			}
		});
		
		panelDestinationDir.add(btnBaseDir);
		txtDestinationDirectory.setBounds(319, 0, 400, 37);
		
		panelDestinationDir.add(txtDestinationDirectory);
		txtcBaseOfDestination.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtcBaseOfDestination);
		txtcBaseOfDestination.setEditable(false);
		txtcBaseOfDestination.setText("Base of Destination Directory");
		txtcBaseOfDestination.setFont(GLSeq2_Main_Application.TEXT_FONT);
		txtcBaseOfDestination.setBounds(0, 0, 261, 37);
		
		panelDestinationDir.add(txtcBaseOfDestination);
		btnPreProcFiles.setIcon(new ImageIcon(Data_And_Library.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
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
				}
			}
		});
		
		panelPreProcessedFiles.add(btnPreProcFiles);
		txtPreProcessedFiles.setBounds(319, 0, 400, 37);
		
		panelPreProcessedFiles.add(txtPreProcessedFiles);
		txtcDirectoryContainingPreprocessed.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtcDirectoryContainingPreprocessed);
		txtcDirectoryContainingPreprocessed.setEditable(false);
		txtcDirectoryContainingPreprocessed.setText("Directory Containing Pre-Processed Files");
		txtcDirectoryContainingPreprocessed.setFont(new Font("Courier New", Font.PLAIN, 11));
		txtcDirectoryContainingPreprocessed.setBounds(0, 0, 261, 37);
		
		panelPreProcessedFiles.add(txtcDirectoryContainingPreprocessed);
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
