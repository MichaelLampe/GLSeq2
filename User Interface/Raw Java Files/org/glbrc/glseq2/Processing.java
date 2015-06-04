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

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JSpinner;

import java.awt.Color;
import javax.swing.SpinnerNumberModel;

public class Processing extends JDialog {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private final JPanel contentPanel = new JPanel();
	private final JSpinner spinMinTrimRead = new JSpinner();
	private final JTextArea txtArtificialName = new JTextArea();
	private final JSpinner spinTrimHead = new JSpinner();
	private final JButton btnTrimReads = new JButton(ButtonEnums.OptionButton.TRIMMING.value);
	private final JPanel panelArtificialSequence = new JPanel();
	private final JTextPane txtcNameOfFasta = new JTextPane();
	private final JPanel panelTrimHead = new JPanel();
	private final JTextPane txtcTrimHead = new JTextPane();
	private final JPanel panelMinimumTrimRead = new JPanel();
	private final JTextPane txtcMinimumTrimRead = new JTextPane();
	private final JTextPane txtchTrimmingAndProcessing = new JTextPane();
	private final JPanel buttonPane = new JPanel();
	private final JButton okButton = new JButton("Apply and Close");
	private final JButton cancelButton = new JButton("Cancel");
	
	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		try {
			Processing dialog = new Processing();
			dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
			dialog.setVisible(true);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Create the dialog.
	 */
	public Processing() {
		initGUI();
		// Assign some variables
		txtArtificialName.setText(GLSeq2_Main_Application.att.getArtificialFASTA());
		spinMinTrimRead.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getMinTrim()));
		spinTrimHead.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getTrimhead()));
		if (GLSeq2_Main_Application.att.getReadTrim().equals("FALSE")){
			btnTrimReads.setText(ButtonEnums.OptionButton.NO_TRIMMING.value);
		}
		
	}
	private void initGUI() {
		setBackground(Color.LIGHT_GRAY);
		setResizable(false);
		setBounds(100, 100, 534, 387);
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBackground(Color.LIGHT_GRAY);
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		contentPanel.setLayout(null);
		{
			panelArtificialSequence.setForeground(Color.DARK_GRAY);
			panelArtificialSequence.setBackground(Color.LIGHT_GRAY);
			panelArtificialSequence.setLayout(null);
			panelArtificialSequence.setBounds(10, 131, 496, 59);
			contentPanel.add(panelArtificialSequence);
			{
				txtArtificialName.setBounds(156, 11, 330, 37);
				panelArtificialSequence.add(txtArtificialName);
			}
			{
				txtcNameOfFasta.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcNameOfFasta);
				txtcNameOfFasta.setText("Name of FASTA File with Artificial Sequences");
				txtcNameOfFasta.setFont(GLSeq2_Main_Application.TEXT_FONT);
				txtcNameOfFasta.setEditable(false);
				txtcNameOfFasta.setBounds(10, 11, 136, 37);
				panelArtificialSequence.add(txtcNameOfFasta);
			}
		}
		{
			panelTrimHead.setForeground(Color.DARK_GRAY);
			panelTrimHead.setBackground(Color.LIGHT_GRAY);
			panelTrimHead.setLayout(null);
			panelTrimHead.setBounds(10, 201, 496, 40);
			contentPanel.add(panelTrimHead);
			{
				spinTrimHead.setModel(new SpinnerNumberModel(new Integer(12), new Integer(0), null, new Integer(1)));
				spinTrimHead.setBounds(116, 0, 98, 40);
				panelTrimHead.add(spinTrimHead);
			}
			{
				txtcTrimHead.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcTrimHead);
				txtcTrimHead.setText("Trim Head");
				txtcTrimHead.setEditable(false);
				txtcTrimHead.setBounds(0, 0, 106, 40);
				panelTrimHead.add(txtcTrimHead);
			}
		}
		{
			panelMinimumTrimRead.setForeground(Color.DARK_GRAY);
			panelMinimumTrimRead.setBackground(Color.LIGHT_GRAY);
			panelMinimumTrimRead.setLayout(null);
			panelMinimumTrimRead.setBounds(10, 252, 496, 40);
			contentPanel.add(panelMinimumTrimRead);
			{
				spinMinTrimRead.setModel(new SpinnerNumberModel(new Integer(36), new Integer(0), null, new Integer(1)));
				spinMinTrimRead.setBounds(116, 0, 98, 40);
				panelMinimumTrimRead.add(spinMinTrimRead);
			}
			{
				txtcMinimumTrimRead.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtcMinimumTrimRead);
				txtcMinimumTrimRead.setText("Minimum Trim Read");
				txtcMinimumTrimRead.setEditable(false);
				txtcMinimumTrimRead.setBounds(0, 0, 106, 40);
				panelMinimumTrimRead.add(txtcMinimumTrimRead);
			}
		}
		{
			btnTrimReads.setForeground(Color.DARK_GRAY);
			btnTrimReads.setFont(new Font("Arial", Font.PLAIN, 20));
			btnTrimReads.setBounds(10, 51, 496, 69);
			btnTrimReads.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent arg0) {
					if (btnTrimReads.getText().equals(ButtonEnums.OptionButton.TRIMMING.value)){
						btnTrimReads.setText(ButtonEnums.OptionButton.NO_TRIMMING.value);
					}
					else{
						btnTrimReads.setText(ButtonEnums.OptionButton.TRIMMING.value);
					}
				}
			});
			contentPanel.add(btnTrimReads);
		}
		{
			nimbusFix(Color.LIGHT_GRAY,txtchTrimmingAndProcessing);
			txtchTrimmingAndProcessing.setForeground(Color.DARK_GRAY);
			txtchTrimmingAndProcessing.setFont(GLSeq2_Main_Application.HEADER_FONT);
			txtchTrimmingAndProcessing.setText("Trimming and Processing Options");
			txtchTrimmingAndProcessing.setEditable(false);
			txtchTrimmingAndProcessing.setBounds(10, 7, 496, 33);
			contentPanel.add(txtchTrimmingAndProcessing);
		}
		{
			buttonPane.setBackground(Color.GRAY);
			buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			{
				okButton.setActionCommand("OK");
				okButton.addActionListener(new ActionListener(){
					public void actionPerformed(ActionEvent arg0){
						GLSeq2_Main_Application.att.setArtificialFASTA(txtArtificialName.getText());
						GLSeq2_Main_Application.att.setMinTrim(String.valueOf(spinMinTrimRead.getValue()));
						GLSeq2_Main_Application.att.setTrimhead(String.valueOf(spinTrimHead.getValue()));
						if (btnTrimReads.getText().equals(ButtonEnums.OptionButton.NO_TRIMMING.value)){
							GLSeq2_Main_Application.att.setReadTrim("FALSE");
						} else{
							GLSeq2_Main_Application.att.setReadTrim("TRUE");
						}
						dispose();
					}
				});
				buttonPane.add(okButton);
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
