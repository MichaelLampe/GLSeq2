package org.glbrc.glseq2;

import java.awt.BorderLayout;
import java.awt.FlowLayout;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.UIDefaults;
import javax.swing.UIManager;
import javax.swing.UIManager.LookAndFeelInfo;
import javax.swing.border.EmptyBorder;
import javax.swing.JTextArea;
import javax.swing.JTextPane;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JSpinner;

import java.awt.Color;
import javax.swing.SpinnerNumberModel;

public class Reference extends JDialog {

	private final JPanel contentPanel = new JPanel();
	private final JTextArea referenceGenome = new JTextArea();
	private final JTextArea refFastaFile = new JTextArea();
	private final JTextArea refGFF = new JTextArea();
	private final JTextArea refFeatureID = new JTextArea();
	private final JSpinner gtfFileColumns = new JSpinner();

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
		referenceGenome.setText(GLSeq2_Main_Application.att.getRGenome());
		refFastaFile.setText(GLSeq2_Main_Application.att.getFASTAname());
		refGFF.setText(GLSeq2_Main_Application.att.getGFFname());
		refFeatureID.setText(GLSeq2_Main_Application.att.getIdAttr());
		gtfFileColumns.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getFeatureColumn()));
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
			JPanel refGenome = new JPanel();
			refGenome.setForeground(Color.DARK_GRAY);
			refGenome.setBackground(Color.LIGHT_GRAY);
			refGenome.setLayout(null);
			refGenome.setBounds(10, 55, 496, 59);
			contentPanel.add(refGenome);
			{
				referenceGenome.setBounds(135, 11, 355, 37);
				refGenome.add(referenceGenome);
			}
			{
				JTextPane txtpnReferenceGenome = new JTextPane();
				txtpnReferenceGenome.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnReferenceGenome);				
				txtpnReferenceGenome.setText("Reference Genome");
				txtpnReferenceGenome.setFont(new Font("Arial", Font.PLAIN, 11));
				txtpnReferenceGenome.setEditable(false);
				txtpnReferenceGenome.setBounds(10, 11, 122, 37);
				refGenome.add(txtpnReferenceGenome);
			}
		}
		{
			JPanel panel = new JPanel();
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setLayout(null);
			panel.setBounds(10, 125, 496, 59);
			contentPanel.add(panel);
			{
				refFastaFile.setBounds(135, 11, 355, 37);
				panel.add(refFastaFile);
			}
			{
				JTextPane txtpnNameOfReference = new JTextPane();
				txtpnNameOfReference.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnNameOfReference);
				txtpnNameOfReference.setText("Name of Reference FASTA File");
				txtpnNameOfReference.setFont(new Font("Arial", Font.PLAIN, 11));
				txtpnNameOfReference.setEditable(false);
				txtpnNameOfReference.setBounds(10, 11, 121, 37);
				panel.add(txtpnNameOfReference);
			}
		}
		{
			JPanel panel = new JPanel();
			panel.setForeground(Color.DARK_GRAY);
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setLayout(null);
			panel.setBounds(10, 195, 496, 59);
			contentPanel.add(panel);
			{
				refGFF.setBounds(135, 11, 355, 37);
				panel.add(refGFF);
			}
			{
				JTextPane txtpnNameOfReference_1 = new JTextPane();
				txtpnNameOfReference_1.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnNameOfReference_1);
				txtpnNameOfReference_1.setText("Name of Reference Genomic Features File");
				txtpnNameOfReference_1.setFont(new Font("Arial", Font.PLAIN, 11));
				txtpnNameOfReference_1.setEditable(false);
				txtpnNameOfReference_1.setBounds(10, 11, 121, 37);
				panel.add(txtpnNameOfReference_1);
			}
		}
		{
			JPanel panel = new JPanel();
			panel.setForeground(Color.DARK_GRAY);
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setLayout(null);
			panel.setBounds(10, 265, 496, 59);
			contentPanel.add(panel);
			{
				refFeatureID.setBounds(135, 11, 355, 37);
				panel.add(refFeatureID);
			}
			{
				JTextPane txtpnNameOfGff = new JTextPane();
				txtpnNameOfGff.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnNameOfGff);
				txtpnNameOfGff.setText("Reference Feature ID");
				txtpnNameOfGff.setFont(new Font("Arial", Font.PLAIN, 11));
				txtpnNameOfGff.setEditable(false);
				txtpnNameOfGff.setBounds(10, 11, 121, 37);
				panel.add(txtpnNameOfGff);
			}
		}
		{
			JTextPane txtpnReferenceFileOptions = new JTextPane();
			txtpnReferenceFileOptions.setFont(new Font("Arial", Font.PLAIN, 20));
			txtpnReferenceFileOptions.setForeground(Color.DARK_GRAY);
			nimbusFix(Color.LIGHT_GRAY,txtpnReferenceFileOptions);
			txtpnReferenceFileOptions.setEditable(false);
			txtpnReferenceFileOptions.setText("Reference File Options");
			txtpnReferenceFileOptions.setBounds(10, 11, 496, 33);
			contentPanel.add(txtpnReferenceFileOptions);
		}
		{
			JTextPane txtpnNumberOfColumns = new JTextPane();
			txtpnNumberOfColumns.setBounds(10, 336, 180, 40);
			contentPanel.add(txtpnNumberOfColumns);
			txtpnNumberOfColumns.setForeground(Color.DARK_GRAY);
			nimbusFix(Color.LIGHT_GRAY,txtpnNumberOfColumns);
			txtpnNumberOfColumns.setText("Number of Columns in the GTF File");
			txtpnNumberOfColumns.setEditable(false);
		}
		{
			gtfFileColumns.setModel(new SpinnerNumberModel(new Integer(9), new Integer(0), null, new Integer(1)));
			gtfFileColumns.setBounds(192, 336, 59, 40);
			contentPanel.add(gtfFileColumns);
		}
		{
			JTextArea txtrOk = new JTextArea();
			txtrOk.setText("OK");
			txtrOk.setBounds(307, 11, 4, 22);
			contentPanel.add(txtrOk);
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
				okButton.addActionListener(new ActionListener(){
					public void actionPerformed(ActionEvent e) {
						GLSeq2_Main_Application.att.setRGenome(referenceGenome.getText());
						GLSeq2_Main_Application.att.setFASTAname(refFastaFile.getText());
						GLSeq2_Main_Application.att.setGFFname(refGFF.getText());
						GLSeq2_Main_Application.att.setIdAttr(refFeatureID.getText());
						GLSeq2_Main_Application.att.setFeatureColumn(String.valueOf(gtfFileColumns.getValue()));
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
