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

import javax.swing.JTextPane;
import javax.swing.text.SimpleAttributeSet;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyledDocument;

public class GLSeq2_Main_Application {

	private JFrame frame;
	public static Attributes att;
	public static RunOptions run;
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
					// Declare new attribute and run
					att = new Attributes();
					run = new RunOptions();
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
		initialize();
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
		
		JButton btnDataAndLibrary = new JButton("Data Sources");
		btnDataAndLibrary.setBounds(33, 54, 300, 70);
		attributeFileContainer.add(btnDataAndLibrary);
		btnDataAndLibrary.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				Data_And_Library data_and_library = new Data_And_Library();
				data_and_library.setVisible(true);
			}
		});
		
		JButton btnPipeline = new JButton("Pipeline");
		btnPipeline.setBounds(33, 135, 300, 70);
		attributeFileContainer.add(btnPipeline);
		btnPipeline.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				Script_Running_Options script_running = new Script_Running_Options();
				script_running.setVisible(true);
			}
		});
		
		JButton btnReference = new JButton("Reference");
		btnReference.setBounds(33, 216, 300, 70);
		attributeFileContainer.add(btnReference);
		btnReference.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				Reference reference = new Reference();
				reference.setVisible(true);
			}
		});
		
		JButton btnProcessing = new JButton("Trimming and Processing");
		btnProcessing.setBounds(33, 297, 300, 70);
		attributeFileContainer.add(btnProcessing);
		btnProcessing.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				Processing processing = new Processing();
				processing.setVisible(true);
			}
		});
		
		JButton btnEnvironment = new JButton("Environment");
		btnEnvironment.setBounds(33, 378, 300, 70);
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
		
		JButton btnGenerateAttributeFile = new JButton("Generate Attribute File");
		btnGenerateAttributeFile.setBounds(33, 491, 300, 70);
		attributeFileContainer.add(btnGenerateAttributeFile);
		
		JPanel panel = new JPanel();
		panel.setBackground(Color.GRAY);
		panel.setLayout(null);
		panel.setBounds(33, 455, 298, 33);
		attributeFileContainer.add(panel);
		
		JTextArea textArea = new JTextArea();
		textArea.setFont(new Font("Arial", Font.PLAIN, 13));
		textArea.setBounds(92, 7, 206, 19);
		panel.add(textArea);
		
		JTextPane txtpnUniqueRunName = new JTextPane();
		txtpnUniqueRunName.setForeground(Color.WHITE);
		txtpnUniqueRunName.setBackground(Color.GRAY);
		txtpnUniqueRunName.setText("Unique Run Name");
		txtpnUniqueRunName.setFont(new Font("Arial", Font.PLAIN, 11));
		txtpnUniqueRunName.setEditable(false);
		txtpnUniqueRunName.setBounds(0, 7, 91, 19);
		panel.add(txtpnUniqueRunName);
		
		JPanel runOptionsContainer = new JPanel();
		runOptionsContainer.setBackground(Color.GRAY);
		runOptionsContainer.setBounds(367, 0, 367, 572);
		frame.getContentPane().add(runOptionsContainer);
		runOptionsContainer.setLayout(null);
		
		JButton btnUpdatingFromDatabase = new JButton("Updating from Database");
		btnUpdatingFromDatabase.setBounds(33, 54, 300, 70);
		runOptionsContainer.add(btnUpdatingFromDatabase);
		
		JButton btnPreprocessingData = new JButton("Pre-Processing Data");
		btnPreprocessingData.setBounds(33, 135, 300, 70);
		runOptionsContainer.add(btnPreprocessingData);
		
		JButton btnAligning = new JButton("Aligning");
		btnAligning.setBounds(33, 216, 300, 70);
		runOptionsContainer.add(btnAligning);
		
		JButton btnCounting = new JButton("Counting");
		btnCounting.setBounds(33, 297, 300, 70);
		runOptionsContainer.add(btnCounting);
		
		JButton btnCollectingResults = new JButton("Collecting Results");
		btnCollectingResults.setBounds(33, 378, 300, 70);
		runOptionsContainer.add(btnCollectingResults);
		
		JButton btnRunningInThe = new JButton("Running in the Background");
		btnRunningInThe.setBounds(33, 459, 300, 70);
		runOptionsContainer.add(btnRunningInThe);
		
		JPanel runOptionsTitleContainer = new JPanel();
		runOptionsTitleContainer.setBackground(Color.GRAY);
		runOptionsTitleContainer.setBounds(0, 0, 367, 70);
		runOptionsContainer.add(runOptionsTitleContainer);
		runOptionsTitleContainer.setLayout(null);
		
		JTextPane runOptionsTitle = new JTextPane();
		runOptionsTitle.setParagraphAttributes(center, false);
		runOptionsTitle.setFont(new Font("Arial", Font.PLAIN, 26));
		runOptionsTitle.setForeground(Color.WHITE);
		runOptionsTitle.setBackground(Color.GRAY);
		runOptionsTitle.setText("Run Options");
		runOptionsTitle.setEditable(false);
		runOptionsTitle.setBounds(10, 5, 347, 48);
		runOptionsTitleContainer.add(runOptionsTitle);
		
		JPanel runContainer = new JPanel();
		runContainer.setBackground(Color.GRAY);
		runContainer.setBounds(734, 0, 361, 572);
		frame.getContentPane().add(runContainer);
		runContainer.setLayout(null);
		
		JPanel imageContainer = new JPanel();
		imageContainer.setBounds(10, 11, 341, 205);
		runContainer.add(imageContainer);
		
		JPanel updateTextStreamContainer = new JPanel();
		updateTextStreamContainer.setBounds(10, 227, 341, 253);
		runContainer.add(updateTextStreamContainer);
		
		JButton btnRun = new JButton("Run");
		btnRun.setBounds(10, 491, 341, 70);
		runContainer.add(btnRun);
	}
}
