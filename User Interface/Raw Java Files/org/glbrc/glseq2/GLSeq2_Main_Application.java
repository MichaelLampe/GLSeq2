package org.glbrc.glseq2;

import java.awt.EventQueue;

import javax.swing.JFrame;
import javax.swing.UIManager;
import javax.swing.UIManager.LookAndFeelInfo;
import javax.swing.JPanel;

import java.awt.Color;

import javax.swing.JMenuBar;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JPopupMenu;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.JButton;
import javax.swing.JTextArea;

import java.awt.Font;
import java.io.IOException;

import javax.swing.JTextPane;
import javax.swing.text.SimpleAttributeSet;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyledDocument;
import javax.swing.JLabel;
import javax.swing.ImageIcon;
import javax.swing.JScrollPane;
import javax.swing.ScrollPaneConstants;

public class GLSeq2_Main_Application {

	public static final JTextPane updates = new JTextPane();;
	private JFrame frame;
	public static Attributes att;
	public static RunOptions run;
	private ScriptTask script;
	private final JTextArea attributeFile = new JTextArea();
	private final JTextPane txtpnAttributeFileLocation = new JTextPane();
	private final JPanel panel = new JPanel();
	private final JScrollPane scrollPane = new JScrollPane();
	private final JButton btnDataAndLibrary = new JButton("Data Sources");
	private final JButton btnPipeline = new JButton("Pipeline");
	private final JButton btnReference = new JButton("Reference");
	private final JButton btnProcessing = new JButton("Trimming and Processing");
	private final JButton btnEnvironment = new JButton("Environment");
	private final JButton generateButton = new JButton("Generate Attribute File");
	private final JButton databaseButton = new JButton("Updating from Database");
	private final JButton preprocessingButton = new JButton("Pre-Processing Data");
	private final JButton alignButton = new JButton("Aligning");
	private final JButton countingButton = new JButton("Counting");
	private final JButton collectingButton = new JButton("Collecting Results");
	private final JButton backgroundButton = new JButton("Running in the Background");
	private final JButton btnRun = new JButton("Run");
	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					GLSeq2_Main_Application window = new GLSeq2_Main_Application();
					window.frame.setVisible(true);
					try {
						for (LookAndFeelInfo info : UIManager
								.getInstalledLookAndFeels()) {
							// Nimbus looks nice, might want to do some
							// improvements to this later.
							if (info.getClassName().contains("Nimbus")) {
								UIManager.setLookAndFeel(info.getClassName());
								break;
							}
						}
					} catch (Exception e) {
						System.out.println("Nimbus not installed.  Using default theme.");
					}
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the application.
	 */
	public GLSeq2_Main_Application() {
		// Declare new attribute and run
		att = new Attributes();
		run = new RunOptions();
		script = new ScriptTask();
		initialize();
		// Updates
		if (run.getUpdateFromDatabase().equals("update")){
			databaseButton.setText("Updating from Database");
		} else{
			databaseButton.setText("NOT Updating from Database");
		}
		//
		if (run.getProcessedData().equals("dataprep")){
			preprocessingButton.setText("Pre-Processing Data");
		} else{
			preprocessingButton.setText("NOT Pre-Processing Data");
		}
		//
		if (run.getAlignment().equals("alignment")){
			alignButton.setText("Aligning");
		} else{
			alignButton.setText("NOT Aligning");
		}
		//
		if (run.getCounting().equals("counting")){
			countingButton.setText("Counting");
		} else{
			countingButton.setText("NOT Counting");
		}
		//
		if (run.getCollectResults().equals("collect")){
			collectingButton.setText("Collecting Results");
		} else{
			collectingButton.setText("NOT Collecting Results");
		}
		//
		if (run.getAmpersand().equals("&")){
			backgroundButton.setText("Running in the Background");
		} else{
			backgroundButton.setText("NOT Running in the Background");
		}
		//	
	}
	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		frame = new JFrame();
		frame.setResizable(false);
		frame.setBounds(100, 100, 1101, 600);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setTitle("GLSeq2 UI");
		frame.getContentPane().setLayout(null);
		
		JPanel attributeFileContainer = new JPanel();
		attributeFileContainer.setBackground(Color.GRAY);
		attributeFileContainer.setBounds(0, 0, 367, 572);
		frame.getContentPane().add(attributeFileContainer);
		attributeFileContainer.setLayout(null);
		
		btnDataAndLibrary.setBounds(33, 40, 300, 70);
		attributeFileContainer.add(btnDataAndLibrary);
		btnDataAndLibrary.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				Data_And_Library data_and_library = new Data_And_Library();
				data_and_library.setVisible(true);
			}
		});
		
		btnPipeline.setBounds(33, 120, 300, 70);
		attributeFileContainer.add(btnPipeline);
		btnPipeline.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				Script_Running_Options script_running = new Script_Running_Options();
				script_running.setVisible(true);
			}
		});
		
		btnReference.setBounds(33, 200, 300, 70);
		attributeFileContainer.add(btnReference);
		btnReference.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				Reference reference = new Reference();
				reference.setVisible(true);
			}
		});
		
		btnProcessing.setBounds(33, 280, 300, 70);
		attributeFileContainer.add(btnProcessing);
		btnProcessing.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				Processing processing = new Processing();
				processing.setVisible(true);
			}
		});
		
		btnEnvironment.setBounds(33, 360, 300, 70);
		attributeFileContainer.add(btnEnvironment);
		btnEnvironment.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				Environment environment = new Environment();
				environment.setVisible(true);
			}
		});
		
		JPanel attributeFileTitleContainer = new JPanel();
		attributeFileTitleContainer.setBackground(Color.GRAY);
		attributeFileTitleContainer.setBounds(0, 0, 367, 70);
		attributeFileContainer.add(attributeFileTitleContainer);
		attributeFileTitleContainer.setLayout(null);
		
		JTextPane attributeFileTitle = new JTextPane();
		//
		StyledDocument doc = attributeFileTitle.getStyledDocument();
		SimpleAttributeSet center = new SimpleAttributeSet();
		StyleConstants.setAlignment(center, StyleConstants.ALIGN_CENTER);
		doc.setParagraphAttributes(0, doc.getLength(), center, false);
		//
		attributeFileTitle.setFont(new Font("Arial", Font.PLAIN, 26));
		attributeFileTitle.setForeground(Color.WHITE);
		attributeFileTitle.setBackground(Color.GRAY);
		attributeFileTitle.setText("Attribute File");
		attributeFileTitle.setEditable(false);
		attributeFileTitle.setBounds(10, 5, 347, 48);
		attributeFileTitleContainer.add(attributeFileTitle);
		generateButton.setBounds(33, 499, 300, 70);
		generateButton.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				try {
					String path = att.writeAttributesFile();
					run.setAttributeFilePath(path);
					attributeFile.setText(path);
				} catch (IOException e1) {
					e1.printStackTrace();
				}
			}
		});
		attributeFileContainer.add(generateButton);
		attributeFile.setBounds(33, 450, 298, 38);
		attributeFileContainer.add(attributeFile);
		attributeFile.setWrapStyleWord(true);
		attributeFile.setLineWrap(true);
		attributeFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtpnAttributeFileLocation.setBounds(135, 430, 91, 19);
		attributeFileContainer.add(txtpnAttributeFileLocation);
		txtpnAttributeFileLocation.setText("Attribute File Path");
		txtpnAttributeFileLocation.setForeground(Color.WHITE);
		txtpnAttributeFileLocation.setFont(new Font("Arial", Font.PLAIN, 11));
		txtpnAttributeFileLocation.setEditable(false);
		txtpnAttributeFileLocation.setBackground(Color.GRAY);
		JPanel runOptionsContainer = new JPanel();
		runOptionsContainer.setBackground(Color.GRAY);
		runOptionsContainer.setBounds(367, 40, 367, 532);
		frame.getContentPane().add(runOptionsContainer);
		runOptionsContainer.setLayout(null);
		databaseButton.setBounds(10, 0, 347, 70);
		databaseButton.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				if (databaseButton.getText().equals("Updating from Database")){
					databaseButton.setText("NOT Updating from Database");
					run.setUpdateFromDatabase("noupdate");
				} else{
					databaseButton.setText("Updating from Database");
					run.setUpdateFromDatabase("update");
				}
			}
		});
		runOptionsContainer.add(databaseButton);
		preprocessingButton.setBounds(10, 81, 347, 70);
		preprocessingButton.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				if (preprocessingButton.getText().equals("Pre-Processing Data")){
					preprocessingButton.setText("NOT Pre-Processing Data");
					run.setProcessedData("dataprep");
				} else{
					preprocessingButton.setText("Pre-Processing Data");
					run.setProcessedData("nodataprep");
				}
			}
		});
		runOptionsContainer.add(preprocessingButton);
		alignButton.setBounds(10, 162, 347, 70);
		alignButton.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				if (alignButton.getText().equals("Aligning")){
					alignButton.setText("NOT Aligning");
					run.setAlignment("noalignment");
				} else{
					alignButton.setText("Aligning");
					run.setAlignment("alignment");
				}
			}
		});
		runOptionsContainer.add(alignButton);
		countingButton.setBounds(10, 243, 347, 70);
		countingButton.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				if (countingButton.getText().equals("Counting")){
					countingButton.setText("NOT Counting");
					run.setCounting("nocounting");
				} else{
					countingButton.setText("Counting");
					run.setCounting("counting");
				}
			}
		});
		runOptionsContainer.add(countingButton);
		collectingButton.setBounds(10, 324, 347, 70);
		collectingButton.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				if (collectingButton.getText().equals("Collecting Results")){
					collectingButton.setText("NOT Collecting Results");
					run.setCollectResults("nocollect");
				} else{
					collectingButton.setText("Collecting Results");
					run.setCollectResults("collect");
				}
			}
		});
		runOptionsContainer.add(collectingButton);
		backgroundButton.setBounds(10, 405, 347, 70);
		backgroundButton.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				if (backgroundButton.getText().equals("Running in the Background")){
					backgroundButton.setText("NOT Running in the Background");
					run.setAmpersand("");
				} else{
					backgroundButton.setText("Running in the Background");
					run.setAmpersand("&");
				}
			}
		});
		runOptionsContainer.add(backgroundButton);
		
		final JTextArea runName = new JTextArea();
		runName.setBounds(111, 495, 246, 19);
		runOptionsContainer.add(runName);
		runName.setFont(new Font("Arial", Font.PLAIN, 13));
		
		JTextPane txtpnUniqueRunName = new JTextPane();
		txtpnUniqueRunName.setBounds(10, 495, 91, 19);
		runOptionsContainer.add(txtpnUniqueRunName);
		txtpnUniqueRunName.setForeground(Color.WHITE);
		txtpnUniqueRunName.setBackground(Color.GRAY);
		txtpnUniqueRunName.setText("Unique Run Name");
		txtpnUniqueRunName.setFont(new Font("Arial", Font.PLAIN, 11));
		txtpnUniqueRunName.setEditable(false);
		
		JPanel runContainer = new JPanel();
		runContainer.setBackground(Color.GRAY);
		runContainer.setBounds(734, 0, 361, 572);
		frame.getContentPane().add(runContainer);
		runContainer.setLayout(null);
		panel.setBackground(Color.LIGHT_GRAY);
		panel.setBounds(10, 40, 341, 440);
		
		runContainer.add(panel);
		panel.setLayout(null);
		scrollPane.setBounds(10, 11, 321, 418);
		
		panel.add(scrollPane);
		updates.setFont(new Font("Arial", Font.PLAIN, 10));
		updates.setEnabled(false);
		updates.setEditable(false);
		
		scrollPane.setViewportView(updates);
		btnRun.setBounds(10, 491, 341, 70);
		btnRun.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				run.setAttributeFilePath(attributeFile.getText());
				run.setRunId(runName.getText());
				updates.setEnabled(true);
				updating("Now running script with arguments: " + String.valueOf(run.returnArgs()));
				script.execute();
			}
		});
		runContainer.add(btnRun);
		
		JTextPane txtpnRunningUpdates = new JTextPane();
		txtpnRunningUpdates.setText("Running Updates");
		StyledDocument docs = txtpnRunningUpdates.getStyledDocument();
		docs.setParagraphAttributes(0, docs.getLength(), center, false);
		txtpnRunningUpdates.setForeground(Color.WHITE);
		txtpnRunningUpdates.setFont(new Font("Arial", Font.PLAIN, 26));
		txtpnRunningUpdates.setEditable(false);
		txtpnRunningUpdates.setBackground(Color.GRAY);
		txtpnRunningUpdates.setBounds(10, 5, 341, 48);
		runContainer.add(txtpnRunningUpdates);
		
		JPanel runOptionsTitleContainer = new JPanel();
		runOptionsTitleContainer.setBounds(367, 0, 379, 213);
		frame.getContentPane().add(runOptionsTitleContainer);
		runOptionsTitleContainer.setBackground(Color.GRAY);
		runOptionsTitleContainer.setLayout(null);
		
		JTextPane runOptionsTitle = new JTextPane();
		runOptionsTitle.setBounds(10, 5, 347, 48);
		runOptionsTitleContainer.add(runOptionsTitle);
		runOptionsTitle.setParagraphAttributes(center, false);
		runOptionsTitle.setFont(new Font("Arial", Font.PLAIN, 26));
		runOptionsTitle.setForeground(Color.WHITE);
		runOptionsTitle.setBackground(Color.GRAY);
		runOptionsTitle.setText("Run Options");
		runOptionsTitle.setEditable(false);
	}
	public static void updating(String the_update){
		String current = updates.getText();
		updates.setText(current + "\n" + the_update);
	}
}
