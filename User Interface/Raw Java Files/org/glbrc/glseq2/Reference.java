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

import javax.swing.JSpinner;

import java.awt.Color;
import javax.swing.SpinnerNumberModel;

public class Reference extends JDialog {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private final JPanel contentPanel = new JPanel();
	private final JTextArea txtReferenceGenome = new JTextArea();
	private final JTextArea txtRefFastaFile = new JTextArea();
	private final JTextArea txtRefGFF = new JTextArea();
	private final JTextArea txtRefFeatureID = new JTextArea();
	private final JSpinner spinGtfFileColumns = new JSpinner();
	private final JPanel panelReferenceGenome = new JPanel();
	private final JTextPane txtcReferenceGenome = new JTextPane();
	private final JPanel panelReferenceFastaFile = new JPanel();
	private final JTextPane txtcNameOfReference = new JTextPane();
	private final JPanel panelGenomicFeaturesFile = new JPanel();
	private final JTextPane txtcRefGenomicFeatures = new JTextPane();
	private final JPanel panelRefFeatureId = new JPanel();
	private final JTextPane txtcNameOfGff = new JTextPane();
	private final JTextPane txtchReferenceFileOptions = new JTextPane();
	private final JTextPane txtcNumberOfColumns = new JTextPane();
	private final JTextArea txtrOk = new JTextArea();
	private final JPanel buttonPane = new JPanel();
	private final JButton okButton = new JButton("Apply and Close");
	private final JButton cancelButton = new JButton("Cancel");
	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		try {
			Reference dialog = new Reference();
			dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
			dialog.setVisible(true);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Create the dialog.
	 */
	public Reference() {
		initGUI();
		//
		txtReferenceGenome.setText(GLSeq2_Main_Application.att.getRGenome());
		txtRefFastaFile.setText(GLSeq2_Main_Application.att.getFASTAname());
		txtRefGFF.setText(GLSeq2_Main_Application.att.getGFFname());
		txtRefFeatureID.setText(GLSeq2_Main_Application.att.getIdAttr());
		spinGtfFileColumns.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getFeatureColumn()));
	}
	private void initGUI() {
		setResizable(false);
		setBounds(100, 100, 536, 463);
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBackground(Color.LIGHT_GRAY);
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		contentPanel.setLayout(null);
		{
			panelReferenceGenome.setForeground(Color.DARK_GRAY);
			panelReferenceGenome.setBackground(Color.LIGHT_GRAY);
			panelReferenceGenome.setLayout(null);
			panelReferenceGenome.setBounds(10, 55, 496, 59);
			contentPanel.add(panelReferenceGenome);
			{
				txtReferenceGenome.setBounds(135, 11, 355, 37);
				panelReferenceGenome.add(txtReferenceGenome);
			}
			{
				txtcReferenceGenome.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcReferenceGenome);				
				txtcReferenceGenome.setText("Reference Genome");
				txtcReferenceGenome.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcReferenceGenome.setEditable(false);
				txtcReferenceGenome.setBounds(10, 11, 122, 37);
				panelReferenceGenome.add(txtcReferenceGenome);
			}
		}
		{
			panelReferenceFastaFile.setBackground(Color.LIGHT_GRAY);
			panelReferenceFastaFile.setLayout(null);
			panelReferenceFastaFile.setBounds(10, 125, 496, 59);
			contentPanel.add(panelReferenceFastaFile);
			{
				txtRefFastaFile.setBounds(135, 11, 355, 37);
				panelReferenceFastaFile.add(txtRefFastaFile);
			}
			{
				txtcNameOfReference.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcNameOfReference);
				txtcNameOfReference.setText("Name of Reference FASTA File");
				txtcNameOfReference.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcNameOfReference.setEditable(false);
				txtcNameOfReference.setBounds(10, 11, 121, 37);
				panelReferenceFastaFile.add(txtcNameOfReference);
			}
		}
		{
			panelGenomicFeaturesFile.setForeground(Color.DARK_GRAY);
			panelGenomicFeaturesFile.setBackground(Color.LIGHT_GRAY);
			panelGenomicFeaturesFile.setLayout(null);
			panelGenomicFeaturesFile.setBounds(10, 195, 496, 59);
			contentPanel.add(panelGenomicFeaturesFile);
			{
				txtRefGFF.setBounds(135, 11, 355, 37);
				panelGenomicFeaturesFile.add(txtRefGFF);
			}
			{
				txtcRefGenomicFeatures.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcRefGenomicFeatures);
				txtcRefGenomicFeatures.setText("Name of Reference Genomic Features File");
				txtcRefGenomicFeatures.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcRefGenomicFeatures.setEditable(false);
				txtcRefGenomicFeatures.setBounds(10, 11, 121, 37);
				panelGenomicFeaturesFile.add(txtcRefGenomicFeatures);
			}
		}
		{
			panelRefFeatureId.setForeground(Color.DARK_GRAY);
			panelRefFeatureId.setBackground(Color.LIGHT_GRAY);
			panelRefFeatureId.setLayout(null);
			panelRefFeatureId.setBounds(10, 265, 496, 59);
			contentPanel.add(panelRefFeatureId);
			{
				txtRefFeatureID.setBounds(135, 11, 355, 37);
				panelRefFeatureId.add(txtRefFeatureID);
			}
			{
				txtcNameOfGff.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcNameOfGff);
				txtcNameOfGff.setText("Reference Feature ID");
				txtcNameOfGff.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcNameOfGff.setEditable(false);
				txtcNameOfGff.setBounds(10, 11, 121, 37);
				panelRefFeatureId.add(txtcNameOfGff);
			}
		}
		{
			txtchReferenceFileOptions.setFont(GLSeq2_Main_Application.HEADER_FONT);
			txtchReferenceFileOptions.setForeground(Color.DARK_GRAY);
			nimbusFix(Color.LIGHT_GRAY,txtchReferenceFileOptions);
			txtchReferenceFileOptions.setEditable(false);
			txtchReferenceFileOptions.setText("Reference File Options");
			txtchReferenceFileOptions.setBounds(10, 11, 496, 33);
			contentPanel.add(txtchReferenceFileOptions);
		}
		{
			txtcNumberOfColumns.setBounds(10, 336, 180, 40);
			contentPanel.add(txtcNumberOfColumns);
			txtcNumberOfColumns.setForeground(Color.DARK_GRAY);
			nimbusFix(Color.LIGHT_GRAY,txtcNumberOfColumns);
			txtcNumberOfColumns.setText("Number of Columns in the GTF File");
			txtcNumberOfColumns.setEditable(false);
		}
		{
			spinGtfFileColumns.setModel(new SpinnerNumberModel(new Integer(9), new Integer(0), null, new Integer(1)));
			spinGtfFileColumns.setBounds(192, 336, 59, 40);
			contentPanel.add(spinGtfFileColumns);
		}
		{
			txtrOk.setText("OK");
			txtrOk.setBounds(307, 11, 4, 22);
			contentPanel.add(txtrOk);
		}
		{
			buttonPane.setBackground(Color.GRAY);
			buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			{
				okButton.setActionCommand("OK");
				buttonPane.add(okButton);
				okButton.addActionListener(new ActionListener(){
					public void actionPerformed(ActionEvent e) {
						GLSeq2_Main_Application.att.setRGenome(txtReferenceGenome.getText());
						GLSeq2_Main_Application.att.setFASTAname(txtRefFastaFile.getText());
						GLSeq2_Main_Application.att.setGFFname(txtRefGFF.getText());
						GLSeq2_Main_Application.att.setIdAttr(txtRefFeatureID.getText());
						GLSeq2_Main_Application.att.setFeatureColumn(String.valueOf(spinGtfFileColumns.getValue()));
						dispose();
					}
				});
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
	void nimbusFix(Color background,JTextPane pane){
		  UIDefaults defaults = new UIDefaults();
		  defaults.put("TextPane[Enabled].backgroundPainter", background);
		  pane.putClientProperty("Nimbus.Overrides", defaults);
		  pane.putClientProperty("Nimbus.Overrides.InheritDefaults", true);
		  pane.setBackground(background);
	}
}
