package org.glbrc.glseq2;

import java.awt.BorderLayout;
import java.awt.FlowLayout;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.UIDefaults;
import javax.swing.border.EmptyBorder;
import javax.swing.JTextArea;
import javax.swing.JTextPane;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.Color;
import java.awt.Font;
import java.io.File;

public class Environment extends JDialog {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private final JPanel contentPanel = new JPanel();


	private final JTextPane txtchEnvironmentOptions = new JTextPane();
	/*/
	 * GLSEQ Script
	 */
	private final JTextArea txtGlSeqDirectory = new JTextArea();
	private final JTextPane txtcGlSeqDirectory = new JTextPane();
	private final JButton btnGlSeqDirectory = new JButton("");
	/*
	 *  Picard tools menu
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
	private final JButton btnBWA = new JButton("");
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
	private final JTextPane txtc_Cushaw = new JTextPane();
	/*
	 * Cushaw index menu
	 */
	private final JPanel panelCushaw_Index = new JPanel();
	private final JTextPane txtcCushaw_Index = new JTextPane();
	private final JTextArea txtCushawIndex = new JTextArea();
	private final JButton btnCushawIndex = new JButton("");
	private final JButton btnCushaw = new JButton("");
	/*
	 * Cushaw GPU menu
	 */
	private final JPanel panelCushaw_GPU = new JPanel();
	private final JTextPane txtcCushaw_GPU = new JTextPane();
	private final JTextArea txtCushaw_GPU = new JTextArea();
	private final JButton btnCushawGPU = new JButton("");
	/*
	 * Movement options
	 */
	private final JPanel panelButton = new JPanel();
	private final JButton okButton = new JButton("Apply and Close");
	private final JButton cancelButton = new JButton("Cancel");
	private final JPanel panel = new JPanel();
	
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
		initGUI();
		// Initialize some variables
		txtTrimmomatic.setText(GLSeq2_Main_Application.att.getTrimPath());
		btnTrimmomatic.setIcon(new ImageIcon(Environment.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
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
				}
			}
		});
		
		panelTrimmomatic.add(btnTrimmomatic);
		txtPicard.setText(GLSeq2_Main_Application.att.getPicardToolsPath());
		btnPicard.setIcon(new ImageIcon(Environment.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
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
				}
			}
		});
		
		panelPicard.add(btnPicard);
		txtFastqc.setText(GLSeq2_Main_Application.att.getFastqcPath());
		btnFastqc.setIcon(new ImageIcon(Environment.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
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
				}
			}
		});
		
		panelFastqc.add(btnFastqc);
		txtBwa.setText(GLSeq2_Main_Application.att.getBwaPath());
		btnBWA.setIcon(new ImageIcon(Environment.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
		btnBWA.setBounds(191, 12, 38, 37);
		btnBWA.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				final JFileChooser chooser = new JFileChooser();
				chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				chooser.showOpenDialog(txtBwa);
				try {
					File file = chooser.getSelectedFile();
					txtBwa.setText(file.getAbsolutePath());
				} catch (NullPointerException e) {
				}
			}
		});
		
		panelBwa.add(btnBWA);
		txtBam2Wig.setText(GLSeq2_Main_Application.att.getBam2WigPath());
		btnBam2Wig.setIcon(new ImageIcon(Environment.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
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
				}
			}
		});
		
		panelBam2Wig.add(btnBam2Wig);
		txtCushaw.setText(GLSeq2_Main_Application.att.getCushawPath());
		btnCushaw.setIcon(new ImageIcon(Environment.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
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
				}
			}
		});
		
		panelCushaw.add(btnCushaw);
		txtCushawIndex.setText(GLSeq2_Main_Application.att.getCushawIndexPath());
		btnCushawIndex.setIcon(new ImageIcon(Environment.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
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
				}
			}
		});
		
		panelCushaw_Index.add(btnCushawIndex);
		txtCushaw_GPU.setText(GLSeq2_Main_Application.att.getCushawGpuPath());
		btnCushawGPU.setIcon(new ImageIcon(Environment.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
		btnCushawGPU.setBounds(191, 11, 38, 37);
		btnCushawGPU.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				final JFileChooser chooser = new JFileChooser();
				chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				chooser.showOpenDialog(txtCushaw_GPU);
				try {
					File file = chooser.getSelectedFile();
					txtCushaw_GPU.setText(file.getAbsolutePath());
				} catch (NullPointerException e) {
				}
			}
		});
		
		panelCushaw_GPU.add(btnCushawGPU);
		panel.setLayout(null);
		panel.setForeground(Color.DARK_GRAY);
		panel.setBackground(Color.LIGHT_GRAY);
		panel.setBounds(10, 48, 803, 59);
		
		contentPanel.add(panel);
		txtGlSeqDirectory.setText(GLSeq2_Main_Application.att.getScriptDirectory());
		txtGlSeqDirectory.setBounds(239, 11, 554, 37);
		
		panel.add(txtGlSeqDirectory);
		txtcGlSeqDirectory.setText("Path to GLSeq Scripts");
		txtcGlSeqDirectory.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtcGlSeqDirectory);
		txtcGlSeqDirectory.setFont(new Font("Monospaced", Font.PLAIN, 11));
		txtcGlSeqDirectory.setEditable(false);
		txtcGlSeqDirectory.setBounds(10, 11, 171, 37);
		
		panel.add(txtcGlSeqDirectory);
		btnGlSeqDirectory.setIcon(new ImageIcon(Environment.class.getResource("/com/sun/java/swing/plaf/windows/icons/TreeOpen.gif")));
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
				}
			}
		});
		panel.add(btnGlSeqDirectory);

	}
	private void initGUI() {
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
			panelFastqc.setBounds(10, 236, 803, 59);
			contentPanel.add(panelFastqc);
			{
				txtFastqc.setBounds(239, 11, 554, 37);
				panelFastqc.add(txtFastqc);
			}
			{
				txtcFastqc.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcFastqc);
				txtcFastqc.setText("Path to Fastqc");
				txtcFastqc.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcFastqc.setEditable(false);
				txtcFastqc.setBounds(10, 11, 171, 37);
				panelFastqc.add(txtcFastqc);
			}
		}
		{
			panelPicard.setForeground(Color.DARK_GRAY);
			panelPicard.setBackground(Color.LIGHT_GRAY);
			panelPicard.setLayout(null);
			panelPicard.setBounds(10, 177, 803, 59);
			contentPanel.add(panelPicard);
			{
				txtPicard.setBounds(239, 11, 554, 37);
				panelPicard.add(txtPicard);
			}
			{
				txtcPicard.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcPicard);
				txtcPicard.setText("Path to PicardTools Jar Directory");
				txtcPicard.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcPicard.setEditable(false);
				txtcPicard.setBounds(10, 11, 171, 37);
				panelPicard.add(txtcPicard);
			}
		}
		{
			panelTrimmomatic.setForeground(Color.DARK_GRAY);
			panelTrimmomatic.setBackground(Color.LIGHT_GRAY);
			panelTrimmomatic.setLayout(null);
			panelTrimmomatic.setBounds(10, 118, 803, 59);
			contentPanel.add(panelTrimmomatic);
			{
				txtTrimmomatic.setBounds(239, 11, 554, 37);
				panelTrimmomatic.add(txtTrimmomatic);
			}
			{
				txtcTrimmomatic.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcTrimmomatic);
				txtcTrimmomatic.setText("Path to Trimmomatic");
				txtcTrimmomatic.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcTrimmomatic.setEditable(false);
				txtcTrimmomatic.setBounds(10, 11, 171, 37);
				panelTrimmomatic.add(txtcTrimmomatic);
			}
		}
		{
			panelCushaw.setForeground(Color.DARK_GRAY);
			panelCushaw.setBackground(Color.LIGHT_GRAY);
			panelCushaw.setLayout(null);
			panelCushaw.setBounds(10, 467, 803, 59);
			contentPanel.add(panelCushaw);
			{
				txtCushaw.setBounds(239, 11, 554, 37);
				panelCushaw.add(txtCushaw);
			}
			{
				txtc_Cushaw.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtc_Cushaw);
				txtc_Cushaw.setText("Path to CUSHAW");
				txtc_Cushaw.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtc_Cushaw.setEditable(false);
				txtc_Cushaw.setBounds(10, 11, 171, 37);
				panelCushaw.add(txtc_Cushaw);
			}
		}
		{
			panelCushaw_GPU.setForeground(Color.DARK_GRAY);
			panelCushaw_GPU.setBackground(Color.LIGHT_GRAY);
			panelCushaw_GPU.setLayout(null);
			panelCushaw_GPU.setBounds(10, 526, 803, 59);
			contentPanel.add(panelCushaw_GPU);
			{
				txtCushaw_GPU.setBounds(239, 11, 554, 37);
				panelCushaw_GPU.add(txtCushaw_GPU);
			}
			{
				txtcCushaw_GPU.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcCushaw_GPU);
				txtcCushaw_GPU.setText("Path to CUSHAW-GPU");
				txtcCushaw_GPU.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcCushaw_GPU.setEditable(false);
				txtcCushaw_GPU.setBounds(10, 11, 171, 37);
				panelCushaw_GPU.add(txtcCushaw_GPU);
			}
		}
		{
			panelBam2Wig.setForeground(Color.DARK_GRAY);
			panelBam2Wig.setBackground(Color.LIGHT_GRAY);
			panelBam2Wig.setLayout(null);
			panelBam2Wig.setBounds(10, 352, 803, 59);
			contentPanel.add(panelBam2Wig);
			{
				txtBam2Wig.setBounds(239, 11, 554, 37);
				panelBam2Wig.add(txtBam2Wig);
			}
			{
				txtcBam2Wig.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcBam2Wig);
				txtcBam2Wig.setText("Path to Bam2Wig Shell Script");
				txtcBam2Wig.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcBam2Wig.setEditable(false);
				txtcBam2Wig.setBounds(10, 11, 171, 37);
				panelBam2Wig.add(txtcBam2Wig);
			}
		}
		{
			panelBwa.setForeground(Color.DARK_GRAY);
			panelBwa.setBackground(Color.LIGHT_GRAY);
			panelBwa.setLayout(null);
			panelBwa.setBounds(10, 295, 803, 59);
			contentPanel.add(panelBwa);
			{
				txtBwa.setBounds(239, 11, 554, 37);
				panelBwa.add(txtBwa);
			}
			{
				txtcBwa.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcBwa);
				txtcBwa.setText("Path to BWA");
				txtcBwa.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcBwa.setEditable(false);
				txtcBwa.setBounds(10, 11, 171, 37);
				panelBwa.add(txtcBwa);
			}
		}
		{
			txtchEnvironmentOptions.setFont(GLSeq2_Main_Application.HEADER_FONT);
			txtchEnvironmentOptions.setForeground(Color.DARK_GRAY);
			nimbusFix(Color.LIGHT_GRAY,txtchEnvironmentOptions);
			txtchEnvironmentOptions.setText("Environment Options");
			txtchEnvironmentOptions.setEditable(false);
			txtchEnvironmentOptions.setBounds(10, 4, 496, 33);
			contentPanel.add(txtchEnvironmentOptions);
		}
		{
			panelCushaw_Index.setLayout(null);
			panelCushaw_Index.setForeground(Color.DARK_GRAY);
			panelCushaw_Index.setBackground(Color.LIGHT_GRAY);
			panelCushaw_Index.setBounds(10, 408, 803, 59);
			contentPanel.add(panelCushaw_Index);
			{
				txtCushawIndex.setBounds(239, 11, 554, 37);
				panelCushaw_Index.add(txtCushawIndex);
			}
			{
				txtcCushaw_Index.setText("Path to CUSHAW Index");
				txtcCushaw_Index.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcCushaw_Index);
				txtcCushaw_Index.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcCushaw_Index.setEditable(false);
				txtcCushaw_Index.setBounds(10, 11, 171, 37);
				panelCushaw_Index.add(txtcCushaw_Index);
			}
		}
		{
			panelButton.setBackground(Color.GRAY);
			panelButton.setLayout(new FlowLayout(FlowLayout.RIGHT));
			getContentPane().add(panelButton, BorderLayout.SOUTH);
			{
				okButton.setActionCommand("OK");
				okButton.addActionListener(new ActionListener(){
					public void actionPerformed(ActionEvent arg0) {
						GLSeq2_Main_Application.att.setTrimPath(txtTrimmomatic.getText());
						GLSeq2_Main_Application.att.setPicardToolsPath(txtPicard.getText());
						GLSeq2_Main_Application.att.setFastqcPath(txtFastqc.getText());
						GLSeq2_Main_Application.att.setBwaPath(txtBwa.getText());
						GLSeq2_Main_Application.att.setBam2WigPath(txtBam2Wig.getText());
						GLSeq2_Main_Application.att.setCushawPath(txtCushaw.getText());
						GLSeq2_Main_Application.att.setCushawIndexPath(txtCushawIndex.getText());
						GLSeq2_Main_Application.att.setCushawGpuPath(txtCushaw_GPU.getText());
						GLSeq2_Main_Application.att.setScriptDirectory(txtGlSeqDirectory.getText());
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
	void nimbusFix(Color background,JTextPane pane){
		  UIDefaults defaults = new UIDefaults();
		  defaults.put("TextPane[Enabled].backgroundPainter", background);
		  pane.putClientProperty("Nimbus.Overrides", defaults);
		  pane.putClientProperty("Nimbus.Overrides.InheritDefaults", true);
		  pane.setBackground(background);
	}
}
