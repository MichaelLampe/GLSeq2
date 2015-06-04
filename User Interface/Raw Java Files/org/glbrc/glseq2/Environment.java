package org.glbrc.glseq2;

import java.awt.BorderLayout;
import java.awt.FlowLayout;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.UIDefaults;
import javax.swing.border.EmptyBorder;
import javax.swing.JTextArea;
import javax.swing.JTextPane;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.Color;

public class Environment extends JDialog {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private final JPanel contentPanel = new JPanel();


	private final JTextPane txtchEnvironmentOptions = new JTextPane();
	/*
	 *  Picard tools menu
	 */
	private final JTextArea txtPicard = new JTextArea();
	private final JPanel panelPicard = new JPanel();
	private final JTextPane txtcPicard = new JTextPane();
	/*
	 * Fastqc menu
	 */
	private final JTextArea txtFastqc = new JTextArea();
	private final JPanel panelFastqc = new JPanel();
	private final JTextPane txtcFastqc = new JTextPane();
	/*
	 * Bam2wig menu
	 */
	private final JTextArea txtBam2Wig = new JTextArea();
	private final JPanel panelBam2Wig = new JPanel();
	private final JTextPane txtcBam2Wig = new JTextPane();
	/*
	 * Bwa menu
	 */
	private final JPanel panelBwa = new JPanel();
	private final JTextPane txtcBwa = new JTextPane();
	private final JTextArea txtBwa = new JTextArea();
	/*
	 * Trimmomatic menu
	 */
	private final JPanel panelTrimmomatic = new JPanel();
	private final JTextArea txtTrimmomatic = new JTextArea();
	private final JTextPane txtcTrimmomatic = new JTextPane();
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
	/*
	 * Cushaw GPU menu
	 */
	private final JPanel panelCushaw_GPU = new JPanel();
	private final JTextPane txtcCushaw_GPU = new JTextPane();
	private final JTextArea txtCushaw_GPU = new JTextArea();
	
	/*
	 * Movement options
	 */
	private final JPanel panelButton = new JPanel();
	private final JButton okButton = new JButton("Apply and Close");
	private final JButton cancelButton = new JButton("Cancel");
	
	
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
		txtPicard.setText(GLSeq2_Main_Application.att.getPicardToolsPath());
		txtFastqc.setText(GLSeq2_Main_Application.att.getFastqcPath());
		txtBwa.setText(GLSeq2_Main_Application.att.getBwaPath());
		txtBam2Wig.setText(GLSeq2_Main_Application.att.getBam2WigPath());
		txtCushaw.setText(GLSeq2_Main_Application.att.getCushawPath());
		txtCushawIndex.setText(GLSeq2_Main_Application.att.getCushawIndexPath());
		txtCushaw_GPU.setText(GLSeq2_Main_Application.att.getCushawGpuPath());

	}
	private void initGUI() {
		setBackground(Color.LIGHT_GRAY);
		setBounds(100, 100, 535, 608);
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBackground(Color.LIGHT_GRAY);
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		contentPanel.setLayout(null);
		{
			panelFastqc.setForeground(Color.DARK_GRAY);
			panelFastqc.setBackground(Color.LIGHT_GRAY);
			panelFastqc.setLayout(null);
			panelFastqc.setBounds(10, 166, 496, 59);
			contentPanel.add(panelFastqc);
			{
				txtFastqc.setBounds(156, 11, 330, 37);
				panelFastqc.add(txtFastqc);
			}
			{
				txtcFastqc.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcFastqc);
				txtcFastqc.setText("Path to Fastqc");
				txtcFastqc.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcFastqc.setEditable(false);
				txtcFastqc.setBounds(10, 11, 136, 37);
				panelFastqc.add(txtcFastqc);
			}
		}
		{
			panelPicard.setForeground(Color.DARK_GRAY);
			panelPicard.setBackground(Color.LIGHT_GRAY);
			panelPicard.setLayout(null);
			panelPicard.setBounds(10, 107, 496, 59);
			contentPanel.add(panelPicard);
			{
				txtPicard.setBounds(156, 11, 330, 37);
				panelPicard.add(txtPicard);
			}
			{
				txtcPicard.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcPicard);
				txtcPicard.setText("Path to PicardTools Jar Directory");
				txtcPicard.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcPicard.setEditable(false);
				txtcPicard.setBounds(10, 11, 136, 37);
				panelPicard.add(txtcPicard);
			}
		}
		{
			panelTrimmomatic.setForeground(Color.DARK_GRAY);
			panelTrimmomatic.setBackground(Color.LIGHT_GRAY);
			panelTrimmomatic.setLayout(null);
			panelTrimmomatic.setBounds(10, 48, 496, 59);
			contentPanel.add(panelTrimmomatic);
			{
				txtTrimmomatic.setBounds(156, 11, 330, 37);
				panelTrimmomatic.add(txtTrimmomatic);
			}
			{
				txtcTrimmomatic.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcTrimmomatic);
				txtcTrimmomatic.setText("Path to Trimmomatic");
				txtcTrimmomatic.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcTrimmomatic.setEditable(false);
				txtcTrimmomatic.setBounds(10, 11, 136, 37);
				panelTrimmomatic.add(txtcTrimmomatic);
			}
		}
		{
			panelCushaw.setForeground(Color.DARK_GRAY);
			panelCushaw.setBackground(Color.LIGHT_GRAY);
			panelCushaw.setLayout(null);
			panelCushaw.setBounds(10, 397, 496, 59);
			contentPanel.add(panelCushaw);
			{
				txtCushaw.setBounds(156, 11, 330, 37);
				panelCushaw.add(txtCushaw);
			}
			{
				txtc_Cushaw.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtc_Cushaw);
				txtc_Cushaw.setText("Path to CUSHAW");
				txtc_Cushaw.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtc_Cushaw.setEditable(false);
				txtc_Cushaw.setBounds(10, 11, 136, 37);
				panelCushaw.add(txtc_Cushaw);
			}
		}
		{
			panelCushaw_GPU.setForeground(Color.DARK_GRAY);
			panelCushaw_GPU.setBackground(Color.LIGHT_GRAY);
			panelCushaw_GPU.setLayout(null);
			panelCushaw_GPU.setBounds(10, 456, 496, 59);
			contentPanel.add(panelCushaw_GPU);
			{
				txtCushaw_GPU.setBounds(156, 11, 330, 37);
				panelCushaw_GPU.add(txtCushaw_GPU);
			}
			{
				txtcCushaw_GPU.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcCushaw_GPU);
				txtcCushaw_GPU.setText("Path to CUSHAW-GPU");
				txtcCushaw_GPU.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcCushaw_GPU.setEditable(false);
				txtcCushaw_GPU.setBounds(10, 11, 136, 37);
				panelCushaw_GPU.add(txtcCushaw_GPU);
			}
		}
		{
			panelBam2Wig.setForeground(Color.DARK_GRAY);
			panelBam2Wig.setBackground(Color.LIGHT_GRAY);
			panelBam2Wig.setLayout(null);
			panelBam2Wig.setBounds(10, 282, 496, 59);
			contentPanel.add(panelBam2Wig);
			{
				txtBam2Wig.setBounds(156, 11, 330, 37);
				panelBam2Wig.add(txtBam2Wig);
			}
			{
				txtcBam2Wig.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcBam2Wig);
				txtcBam2Wig.setText("Path to Bam2Wig Shell Script");
				txtcBam2Wig.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcBam2Wig.setEditable(false);
				txtcBam2Wig.setBounds(10, 11, 136, 37);
				panelBam2Wig.add(txtcBam2Wig);
			}
		}
		{
			panelBwa.setForeground(Color.DARK_GRAY);
			panelBwa.setBackground(Color.LIGHT_GRAY);
			panelBwa.setLayout(null);
			panelBwa.setBounds(10, 225, 496, 59);
			contentPanel.add(panelBwa);
			{
				txtBwa.setBounds(156, 11, 330, 37);
				panelBwa.add(txtBwa);
			}
			{
				txtcBwa.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcBwa);
				txtcBwa.setText("Path to BWA");
				txtcBwa.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcBwa.setEditable(false);
				txtcBwa.setBounds(10, 11, 136, 37);
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
			panelCushaw_Index.setBounds(10, 338, 496, 59);
			contentPanel.add(panelCushaw_Index);
			{
				txtCushawIndex.setBounds(156, 11, 330, 37);
				panelCushaw_Index.add(txtCushawIndex);
			}
			{
				txtcCushaw_Index.setText("Path to CUSHAW Index");
				txtcCushaw_Index.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcCushaw_Index);
				txtcCushaw_Index.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcCushaw_Index.setEditable(false);
				txtcCushaw_Index.setBounds(10, 11, 136, 37);
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
