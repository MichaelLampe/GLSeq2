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
import java.awt.Color;

public class Environment extends JDialog {

	private final JPanel contentPanel = new JPanel();
	private final JTextArea trimmomatic = new JTextArea();
	private final JTextArea fastqc = new JTextArea();
	private final JTextArea picard = new JTextArea();
	private final JTextArea cushaw_index = new JTextArea();
	private final JTextArea cushaw = new JTextArea();
	private final JTextArea bam2wig = new JTextArea();
	private final JTextArea bwa = new JTextArea();
	private final JTextArea cushaw_gpu = new JTextArea();

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
		initGUI();
		// Initialize some variables
		trimmomatic.setText(GLSeq2_Main_Application.att.getTrimPath());
		picard.setText(GLSeq2_Main_Application.att.getPicardToolsPath());
		fastqc.setText(GLSeq2_Main_Application.att.getFastqcPath());
		bwa.setText(GLSeq2_Main_Application.att.getBwaPath());
		bam2wig.setText(GLSeq2_Main_Application.att.getBam2WigPath());
		cushaw.setText(GLSeq2_Main_Application.att.getCushawPath());
		cushaw_index.setText(GLSeq2_Main_Application.att.getCushawIndexPath());
		cushaw_gpu.setText(GLSeq2_Main_Application.att.getCushawGpuPath());

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
			JPanel panel = new JPanel();
			panel.setForeground(Color.DARK_GRAY);
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setLayout(null);
			panel.setBounds(10, 166, 496, 59);
			contentPanel.add(panel);
			{
				fastqc.setBounds(156, 11, 330, 37);
				panel.add(fastqc);
			}
			{
				JTextPane txtpnPathToFastqc = new JTextPane();
				txtpnPathToFastqc.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnPathToFastqc);
				txtpnPathToFastqc.setText("Path to Fastqc");
				txtpnPathToFastqc.setFont(new Font("Arial", Font.PLAIN, 11));
				txtpnPathToFastqc.setEditable(false);
				txtpnPathToFastqc.setBounds(10, 11, 136, 37);
				panel.add(txtpnPathToFastqc);
			}
		}
		{
			JPanel panel = new JPanel();
			panel.setForeground(Color.DARK_GRAY);
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setLayout(null);
			panel.setBounds(10, 107, 496, 59);
			contentPanel.add(panel);
			{
				picard.setBounds(156, 11, 330, 37);
				panel.add(picard);
			}
			{
				JTextPane txtpnPathToPicardtools = new JTextPane();
				txtpnPathToPicardtools.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnPathToPicardtools);
				txtpnPathToPicardtools.setText("Path to PicardTools Jar Directory");
				txtpnPathToPicardtools.setFont(new Font("Arial", Font.PLAIN, 11));
				txtpnPathToPicardtools.setEditable(false);
				txtpnPathToPicardtools.setBounds(10, 11, 136, 37);
				panel.add(txtpnPathToPicardtools);
			}
		}
		{
			JPanel panel = new JPanel();
			panel.setForeground(Color.DARK_GRAY);
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setLayout(null);
			panel.setBounds(10, 48, 496, 59);
			contentPanel.add(panel);
			{
				trimmomatic.setBounds(156, 11, 330, 37);
				panel.add(trimmomatic);
			}
			{
				JTextPane txtpnPathToTrimmomatic = new JTextPane();
				txtpnPathToTrimmomatic.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnPathToTrimmomatic);
				txtpnPathToTrimmomatic.setText("Path to Trimmomatic");
				txtpnPathToTrimmomatic.setFont(new Font("Arial", Font.PLAIN, 11));
				txtpnPathToTrimmomatic.setEditable(false);
				txtpnPathToTrimmomatic.setBounds(10, 11, 136, 37);
				panel.add(txtpnPathToTrimmomatic);
			}
		}
		{
			JPanel panel = new JPanel();
			panel.setForeground(Color.DARK_GRAY);
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setLayout(null);
			panel.setBounds(10, 397, 496, 59);
			contentPanel.add(panel);
			{
				cushaw.setBounds(156, 11, 330, 37);
				panel.add(cushaw);
			}
			{
				JTextPane txtpnPathToCushaw = new JTextPane();
				txtpnPathToCushaw.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnPathToCushaw);
				txtpnPathToCushaw.setText("Path to CUSHAW");
				txtpnPathToCushaw.setFont(new Font("Arial", Font.PLAIN, 11));
				txtpnPathToCushaw.setEditable(false);
				txtpnPathToCushaw.setBounds(10, 11, 136, 37);
				panel.add(txtpnPathToCushaw);
			}
		}
		{
			JPanel panel = new JPanel();
			panel.setForeground(Color.DARK_GRAY);
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setLayout(null);
			panel.setBounds(10, 456, 496, 59);
			contentPanel.add(panel);
			{
				cushaw_gpu.setBounds(156, 11, 330, 37);
				panel.add(cushaw_gpu);
			}
			{
				JTextPane txtpnPathToCushawgpu = new JTextPane();
				txtpnPathToCushawgpu.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnPathToCushawgpu);
				txtpnPathToCushawgpu.setText("Path to CUSHAW-GPU");
				txtpnPathToCushawgpu.setFont(new Font("Arial", Font.PLAIN, 11));
				txtpnPathToCushawgpu.setEditable(false);
				txtpnPathToCushawgpu.setBounds(10, 11, 136, 37);
				panel.add(txtpnPathToCushawgpu);
			}
		}
		{
			JPanel panel = new JPanel();
			panel.setForeground(Color.DARK_GRAY);
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setLayout(null);
			panel.setBounds(10, 282, 496, 59);
			contentPanel.add(panel);
			{
				bam2wig.setBounds(156, 11, 330, 37);
				panel.add(bam2wig);
			}
			{
				JTextPane txtpnPathToBamwig = new JTextPane();
				txtpnPathToBamwig.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnPathToBamwig);
				txtpnPathToBamwig.setText("Path to Bam2Wig Shell Script");
				txtpnPathToBamwig.setFont(new Font("Arial", Font.PLAIN, 11));
				txtpnPathToBamwig.setEditable(false);
				txtpnPathToBamwig.setBounds(10, 11, 136, 37);
				panel.add(txtpnPathToBamwig);
			}
		}
		{
			JPanel panel = new JPanel();
			panel.setForeground(Color.DARK_GRAY);
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setLayout(null);
			panel.setBounds(10, 225, 496, 59);
			contentPanel.add(panel);
			{
				bwa.setBounds(156, 11, 330, 37);
				panel.add(bwa);
			}
			{
				JTextPane txtpnPathToBwa = new JTextPane();
				txtpnPathToBwa.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnPathToBwa);
				txtpnPathToBwa.setText("Path to BWA");
				txtpnPathToBwa.setFont(new Font("Arial", Font.PLAIN, 11));
				txtpnPathToBwa.setEditable(false);
				txtpnPathToBwa.setBounds(10, 11, 136, 37);
				panel.add(txtpnPathToBwa);
			}
		}
		{
			JTextPane txtpnEnvironmentOptions = new JTextPane();
			txtpnEnvironmentOptions.setFont(new Font("Arial", Font.PLAIN, 20));
			txtpnEnvironmentOptions.setForeground(Color.DARK_GRAY);
			nimbusFix(Color.LIGHT_GRAY,txtpnEnvironmentOptions);
			txtpnEnvironmentOptions.setText("Environment Options");
			txtpnEnvironmentOptions.setEditable(false);
			txtpnEnvironmentOptions.setBounds(10, 4, 496, 33);
			contentPanel.add(txtpnEnvironmentOptions);
		}
		{
			JPanel panel = new JPanel();
			panel.setLayout(null);
			panel.setForeground(Color.DARK_GRAY);
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setBounds(10, 338, 496, 59);
			contentPanel.add(panel);
			{
				cushaw_index.setBounds(156, 11, 330, 37);
				panel.add(cushaw_index);
			}
			{
				JTextPane txtpnPathToCushaw_1 = new JTextPane();
				txtpnPathToCushaw_1.setText("Path to CUSHAW Index");
				txtpnPathToCushaw_1.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY,txtpnPathToCushaw_1);
				txtpnPathToCushaw_1.setFont(new Font("Arial", Font.PLAIN, 11));
				txtpnPathToCushaw_1.setEditable(false);
				txtpnPathToCushaw_1.setBounds(10, 11, 136, 37);
				panel.add(txtpnPathToCushaw_1);
			}
		}
		{
			JPanel buttonPane = new JPanel();
			buttonPane.setBackground(Color.GRAY);
			buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			{
				JButton okButton = new JButton("Apply and Close");
				okButton.setActionCommand("OK");
				okButton.addActionListener(new ActionListener(){
					public void actionPerformed(ActionEvent arg0) {
						GLSeq2_Main_Application.att.setTrimPath(trimmomatic.getText());
						GLSeq2_Main_Application.att.setPicardToolsPath(picard.getText());
						GLSeq2_Main_Application.att.setFastqcPath(fastqc.getText());
						GLSeq2_Main_Application.att.setBwaPath(bwa.getText());
						GLSeq2_Main_Application.att.setBam2WigPath(bam2wig.getText());
						GLSeq2_Main_Application.att.setCushawPath(cushaw.getText());
						GLSeq2_Main_Application.att.setCushawIndexPath(cushaw_index.getText());
						GLSeq2_Main_Application.att.setCushawGpuPath(cushaw_gpu.getText());
						dispose();
					}
				});	
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
