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

	private final JPanel contentPanel = new JPanel();

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
			JPanel panel = new JPanel();
			panel.setForeground(Color.DARK_GRAY);
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setLayout(null);
			panel.setBounds(10, 131, 496, 59);
			contentPanel.add(panel);
			{
				JTextArea textArea = new JTextArea();
				textArea.setBounds(156, 11, 330, 37);
				panel.add(textArea);
			}
			{
				JTextPane txtpnNameOfFasta = new JTextPane();
				txtpnNameOfFasta.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnNameOfFasta);
				txtpnNameOfFasta.setText("Name of FASTA File with Artificial Sequences");
				txtpnNameOfFasta.setFont(new Font("Arial", Font.PLAIN, 11));
				txtpnNameOfFasta.setEditable(false);
				txtpnNameOfFasta.setBounds(10, 11, 136, 37);
				panel.add(txtpnNameOfFasta);
			}
		}
		{
			JPanel panel = new JPanel();
			panel.setForeground(Color.DARK_GRAY);
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setLayout(null);
			panel.setBounds(10, 201, 496, 40);
			contentPanel.add(panel);
			{
				JSpinner spinner = new JSpinner();
				spinner.setModel(new SpinnerNumberModel(new Integer(12), new Integer(0), null, new Integer(1)));
				spinner.setBounds(116, 0, 98, 40);
				panel.add(spinner);
			}
			{
				JTextPane txtpnTrimHead = new JTextPane();
				txtpnTrimHead.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnTrimHead);
				txtpnTrimHead.setText("Trim Head");
				txtpnTrimHead.setEditable(false);
				txtpnTrimHead.setBounds(0, 0, 106, 40);
				panel.add(txtpnTrimHead);
			}
		}
		{
			JPanel panel = new JPanel();
			panel.setForeground(Color.DARK_GRAY);
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setLayout(null);
			panel.setBounds(10, 252, 496, 40);
			contentPanel.add(panel);
			{
				JSpinner spinner = new JSpinner();
				spinner.setModel(new SpinnerNumberModel(new Integer(36), new Integer(0), null, new Integer(1)));
				spinner.setBounds(116, 0, 98, 40);
				panel.add(spinner);
			}
			{
				JTextPane txtpnMinimumTrimRead = new JTextPane();
				txtpnMinimumTrimRead.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnMinimumTrimRead);
				txtpnMinimumTrimRead.setText("Minimum Trim Read");
				txtpnMinimumTrimRead.setEditable(false);
				txtpnMinimumTrimRead.setBounds(0, 0, 106, 40);
				panel.add(txtpnMinimumTrimRead);
			}
		}
		{
			final JButton btnNewButton = new JButton("Trimming Raw Reads");
			btnNewButton.setForeground(Color.DARK_GRAY);
			btnNewButton.setFont(new Font("Arial", Font.PLAIN, 20));
			btnNewButton.setBounds(10, 51, 496, 69);
			btnNewButton.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent arg0) {
					if (btnNewButton.getText().equals("Trimming Raw Reads")){
						btnNewButton.setText("Not Trimming Raw Reads");
					}
					else{
						btnNewButton.setText("Trimming Raw Reads");
					}
				}
			});
			contentPanel.add(btnNewButton);
		}
		{
			JTextPane txtpnTrimmingAndProcessing = new JTextPane();
			nimbusFix(Color.LIGHT_GRAY,txtpnTrimmingAndProcessing);
			txtpnTrimmingAndProcessing.setForeground(Color.DARK_GRAY);
			txtpnTrimmingAndProcessing.setFont(new Font("Arial", Font.PLAIN, 20));
			txtpnTrimmingAndProcessing.setText("Trimming and Processing Options");
			txtpnTrimmingAndProcessing.setEditable(false);
			txtpnTrimmingAndProcessing.setBounds(10, 7, 496, 33);
			contentPanel.add(txtpnTrimmingAndProcessing);
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
