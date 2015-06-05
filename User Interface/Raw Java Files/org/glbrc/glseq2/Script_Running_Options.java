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
import javax.swing.JComboBox;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JCheckBox;

public class Script_Running_Options extends JDialog {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private final JPanel contentPanel = new JPanel();
	//
	private final JSpinner spinNumberOfCores = new JSpinner();
	private final JComboBox<String> comboAlignmentAlgo = new JComboBox<String>();
	private final JCheckBox checkFeatureCount = new JCheckBox("Feature Counts");
	private final JCheckBox checkHtSeq = new JCheckBox("HTSeq");
	private final JCheckBox checkRsem = new JCheckBox("RSEM");
	private final JTextArea txtProtocolID = new JTextArea();
	private final JSpinner spinParallelExpression = new JSpinner();
	private final JSpinner spinParallelDataPrep = new JSpinner();
	private final JSpinner spinMaxFragLen = new JSpinner();
	private final JSpinner spinMaxBuffer = new JSpinner();
	private final JButton btnExtractCoverage = new JButton();
	private final JButton btnComputeConfIntervals = new JButton();
	private final JButton btnOutputGenomeBam = new JButton();
	private final JPanel panelGlow = new JPanel();
	private final JPanel panelManualHolder = new JPanel();
	private final JPanel panelProtocol = new JPanel();
	private final JTextPane txtcProtocolId = new JTextPane();
	private final JPanel panelAlignment = new JPanel();
	private final JTextPane txtcAlignmentAlgorithm = new JTextPane();
	private final JPanel panelCounting = new JPanel();
	private final JTextPane txtcCountingMethods = new JTextPane();
	private final JTextPane txtchPipelineOptions = new JTextPane();
	private final JPanel panelNumberOfCores = new JPanel();
	private final JTextPane txtcNumberOfCores = new JTextPane();
	private final JPanel panelParallelExpressionComp = new JPanel();
	private final JTextPane txtcParallelComputationStreams = new JTextPane();
	private final JPanel panelParallelDataComp = new JPanel();
	private final JTextPane txtcParallelComputationStreams_1 = new JTextPane();
	private final JPanel panelMaxFragLen = new JPanel();
	private final JTextPane txtcMaximalFragmentLength = new JTextPane();
	private final JPanel panelMaxBufferConfInts = new JPanel();
	private final JTextPane txtcMaximumAuxiliaryBuffer = new JTextPane();
	private final JPanel buttonPane = new JPanel();
	private final JButton okButton = new JButton("Apply and Close");
	private final JButton cancelButton = new JButton("Cancel");
	private final JPanel panel = new JPanel();
	private final JPanel panel_1 = new JPanel();
	private final JPanel panelReferenceGenome = new JPanel();
	private final JTextArea txtReferenceGenome = new JTextArea();
	private final JTextPane txtcReferenceGenome = new JTextPane();
	private final JPanel panelReferenceFasta = new JPanel();
	private final JTextArea txtReferenceFasta = new JTextArea();
	private final JTextPane txtcReferenceFasta = new JTextPane();
	private final JPanel panelReferenceFeatures = new JPanel();
	private final JTextArea txtReferenceFeatures = new JTextArea();
	private final JTextPane txtcReferenceFeatures = new JTextPane();
	private final JPanel panelReferenceID = new JPanel();
	private final JTextArea txtReferenceID = new JTextArea();
	private final JTextPane txtcReferenceID = new JTextPane();
	private final JTextPane txtcColumnGtf = new JTextPane();
	private final JSpinner spinColumnGtf = new JSpinner();
	private final JTextPane txtpnGlowLoginGoing = new JTextPane();
	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		try {
			Script_Running_Options dialog = new Script_Running_Options();
			dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
			dialog.setVisible(true);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Create the dialog.
	 */
	public Script_Running_Options() {
		initGUI();
		spinNumberOfCores.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getCores()));
		spinParallelExpression.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getStreams()));
		spinParallelDataPrep.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getStreamsDataPrep()));
		spinMaxFragLen.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getFragMaxLength()));
		spinMaxBuffer.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getCiMem()));
		txtReferenceGenome.setText(GLSeq2_Main_Application.att.getRGenome());
		System.out.println(GLSeq2_Main_Application.att.getRGenome());
		txtReferenceFasta.setText(GLSeq2_Main_Application.att.getFASTAname());
		txtReferenceFeatures.setText(GLSeq2_Main_Application.att.getGFFname());
		txtReferenceID.setText(GLSeq2_Main_Application.att.getIdAttr());
		spinColumnGtf.setValue(Integer.valueOf(GLSeq2_Main_Application.att.getFeatureColumn()));
		if (!String.valueOf(comboAlignmentAlgo.getSelectedItem()).contains("Bowtie")) {
			checkRsem.setSelected(false);
			checkRsem.setEnabled(false);
		}
		if (String.valueOf(comboAlignmentAlgo.getSelectedItem()).contains("CUSHAW")){
			if (GLSeq2_Main_Application.att.getGpuAccel().equals("TRUE")){
				comboAlignmentAlgo.setSelectedItem("Cushaw-GPU");
			}
		}
		if(GLSeq2_Main_Application.att.getFeatureCounts().contains("FeatureCounts")){
			checkFeatureCount.setSelected(true);
		}
		if(GLSeq2_Main_Application.att.getHTSeq().contains("HTSeq")){
			checkHtSeq.setSelected(true);
		}
		if(GLSeq2_Main_Application.att.getRSEM().contains("RSEM") && String.valueOf(comboAlignmentAlgo.getSelectedItem()).contains("Bowtie")){
			checkRsem.setSelected(true);
		}
		{
			panel.setLayout(null);
			panel.setForeground(Color.DARK_GRAY);
			panel.setBackground(Color.LIGHT_GRAY);
			panel.setBounds(0, 516, 695, 59);
			panelManualHolder.add(panel);
		}
		{
			panelReferenceGenome.setLayout(null);
			panelReferenceGenome.setForeground(Color.DARK_GRAY);
			panelReferenceGenome.setBackground(Color.LIGHT_GRAY);
			panelReferenceGenome.setBounds(10, 0, 983, 40);
			panelManualHolder.add(panelReferenceGenome);
		}
		{
			txtReferenceGenome.setBounds(254, 0, 719, 37);
			panelReferenceGenome.add(txtReferenceGenome);
		}
		{
			txtcReferenceGenome.setText("Reference Genome");
			txtcReferenceGenome.setForeground(Color.DARK_GRAY);
			txtcReferenceGenome.setFont(new Font("Monospaced", Font.PLAIN, 11));
			txtcReferenceGenome.setEditable(false);
			nimbusFix(Color.LIGHT_GRAY,txtcReferenceGenome);	
			txtcReferenceGenome.setBounds(0, 0, 235, 37);
			panelReferenceGenome.add(txtcReferenceGenome);
		}
		{
			panelReferenceFasta.setLayout(null);
			panelReferenceFasta.setBackground(Color.LIGHT_GRAY);
			panelReferenceFasta.setBounds(10, 41, 983, 40);
			panelManualHolder.add(panelReferenceFasta);
		}
		{
			txtReferenceFasta.setBounds(255, 0, 718, 37);
			panelReferenceFasta.add(txtReferenceFasta);
		}
		{
			txtcReferenceFasta.setText("Name of Reference FASTA File");
			txtcReferenceFasta.setForeground(Color.DARK_GRAY);
			txtcReferenceFasta.setFont(new Font("Monospaced", Font.PLAIN, 11));
			nimbusFix(Color.LIGHT_GRAY,txtcReferenceFasta);
			txtcReferenceFasta.setEditable(false);
			txtcReferenceFasta.setBounds(0, 0, 235, 37);
			panelReferenceFasta.add(txtcReferenceFasta);
		}
		{
			panelReferenceFeatures.setLayout(null);
			panelReferenceFeatures.setForeground(Color.DARK_GRAY);
			panelReferenceFeatures.setBackground(Color.LIGHT_GRAY);
			panelReferenceFeatures.setBounds(10, 81, 983, 40);
			panelManualHolder.add(panelReferenceFeatures);
		}
		{
			txtReferenceFeatures.setBounds(255, 0, 718, 37);
			panelReferenceFeatures.add(txtReferenceFeatures);
		}
		{
			txtcReferenceFeatures.setText("Name of Reference Genomic Features File");
			txtcReferenceFeatures.setForeground(Color.DARK_GRAY);
			nimbusFix(Color.LIGHT_GRAY,txtcReferenceFeatures);
			txtcReferenceFeatures.setFont(new Font("Monospaced", Font.PLAIN, 11));
			txtcReferenceFeatures.setEditable(false);
			txtcReferenceFeatures.setBounds(0, 0, 235, 37);
			panelReferenceFeatures.add(txtcReferenceFeatures);
		}
		{
			panelReferenceID.setLayout(null);
			panelReferenceID.setForeground(Color.DARK_GRAY);
			panelReferenceID.setBackground(Color.LIGHT_GRAY);
			panelReferenceID.setBounds(10, 122, 983, 47);
			panelManualHolder.add(panelReferenceID);
		}
		{
			txtReferenceID.setBounds(255, 0, 718, 37);
			panelReferenceID.add(txtReferenceID);
		}
		{
			txtcReferenceID.setText("Reference Feature ID");
			txtcReferenceID.setForeground(Color.DARK_GRAY);
			nimbusFix(Color.LIGHT_GRAY,txtcReferenceID);
			txtcReferenceID.setFont(new Font("Monospaced", Font.PLAIN, 11));
			txtcReferenceID.setEditable(false);
			txtcReferenceID.setBounds(0, 0, 235, 37);
			panelReferenceID.add(txtcReferenceID);
		}
		txtcColumnGtf.setText("Number of Columns in the GTF File");
		txtcColumnGtf.setForeground(Color.DARK_GRAY);
		nimbusFix(Color.LIGHT_GRAY,txtcColumnGtf);
		txtcColumnGtf.setEditable(false);
		txtcColumnGtf.setBounds(550, 379, 368, 40);
		
		panelManualHolder.add(txtcColumnGtf);
		spinColumnGtf.setBounds(929, 379, 53, 40);
		
		panelManualHolder.add(spinColumnGtf);
		{
			panel_1.setLayout(null);
			panel_1.setBackground(Color.LIGHT_GRAY);
			panel_1.setBounds(314, 0, 699, 180);
			contentPanel.add(panel_1);
		}
		txtProtocolID.setText(GLSeq2_Main_Application.run.getProtocolId());
		comboAlignmentAlgo.setSelectedItem(GLSeq2_Main_Application.att.getaAlgor());
		
		if (GLSeq2_Main_Application.att.getStrandExtract().equals("TRUE")){
			btnExtractCoverage.setText(ButtonEnums.OptionButton.EXTRACT.value);
		} else{
			btnExtractCoverage.setText(ButtonEnums.OptionButton.NO_EXTRACT.value);
		}
		
		if (GLSeq2_Main_Application.att.getCompConf().equals("TRUE")){
			btnComputeConfIntervals.setText(ButtonEnums.OptionButton.COMPUTE.value);
		} else{
			btnComputeConfIntervals.setText(ButtonEnums.OptionButton.NO_COMPUTE.value);
		}
		
		if (GLSeq2_Main_Application.att.getGenoBam().equals("TRUE")){
			btnOutputGenomeBam.setText(ButtonEnums.OptionButton.OUTPUT.value);
		} else{
			btnOutputGenomeBam.setText(ButtonEnums.OptionButton.NO_OUTPUT.value);
		}
	}

	private void initGUI() {
		setBounds(100, 100, 1029, 693);
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBackground(Color.LIGHT_GRAY);
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		contentPanel.setLayout(null);
		{
			panelGlow.setBackground(Color.ORANGE);
			panelGlow.setLayout(null);
			panelGlow.setBounds(0, 0, 298, 180);
			contentPanel.add(panelGlow);
		}
		{
			txtpnGlowLoginGoing.setText("Glow Login Going Here");
			txtpnGlowLoginGoing.setBounds(39, 38, 204, 20);
			panelGlow.add(txtpnGlowLoginGoing);
		}
		{
			panelManualHolder.setBackground(Color.LIGHT_GRAY);
			panelManualHolder.setBounds(10, 191, 1003, 525);
			contentPanel.add(panelManualHolder);
			panelManualHolder.setLayout(null);
			{
				panelProtocol.setBounds(295, 11, 394, 59);
				panel_1.add(panelProtocol);
				panelProtocol.setForeground(Color.DARK_GRAY);
				panelProtocol.setBackground(Color.LIGHT_GRAY);
				panelProtocol.setLayout(null);
				{
					txtProtocolID.setBounds(181, 11, 203, 37);
					panelProtocol.add(txtProtocolID);
				}
				{
					txtcProtocolId.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtcProtocolId);
					txtcProtocolId.setText("Protocol ID");
					txtcProtocolId.setFont(GLSeq2_Main_Application.TEXT_FONT);
					txtcProtocolId.setEditable(false);
					txtcProtocolId.setBounds(10, 11, 161, 37);
					panelProtocol.add(txtcProtocolId);
				}
			}
			{
				panelAlignment.setBounds(10, 81, 354, 92);
				panel_1.add(panelAlignment);
				panelAlignment.setForeground(Color.DARK_GRAY);
				panelAlignment.setBackground(Color.LIGHT_GRAY);
				panelAlignment.setLayout(null);
				{
					txtcAlignmentAlgorithm.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtcAlignmentAlgorithm);
					txtcAlignmentAlgorithm.setText("Alignment Algorithm");
					txtcAlignmentAlgorithm.setEditable(false);
					txtcAlignmentAlgorithm.setBounds(10, 0, 334, 20);
					panelAlignment.add(txtcAlignmentAlgorithm);
				}
				comboAlignmentAlgo.setFont(GLSeq2_Main_Application.TEXT_FONT);
				comboAlignmentAlgo.setModel(new DefaultComboBoxModel<String>(
						new String[] { "BWA", "Bowtie", "Bowtie2", "Cushaw",
								"Cushaw-GPU" }));
				comboAlignmentAlgo.setBounds(10, 28, 334, 42);
				comboAlignmentAlgo.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						if (String.valueOf(comboAlignmentAlgo.getSelectedItem())
								.contains("Bowtie")) {
							checkRsem.setEnabled(true);
						} else {
							checkRsem.setSelected(false);
							checkRsem.setEnabled(false);
						}
					}
				});
				panelAlignment.add(comboAlignmentAlgo);
			}
			{
				panelCounting.setBounds(368, 81, 321, 92);
				panel_1.add(panelCounting);
				panelCounting.setForeground(Color.DARK_GRAY);
				panelCounting.setBackground(Color.LIGHT_GRAY);
				panelCounting.setLayout(null);
				{
					txtcCountingMethods.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtcCountingMethods);
					txtcCountingMethods.setText("Counting Methods");
					txtcCountingMethods.setEditable(false);
					txtcCountingMethods.setBounds(10, 0, 301, 20);
					panelCounting.add(txtcCountingMethods);
				}
				checkFeatureCount.setFont(GLSeq2_Main_Application.TEXT_FONT);
				checkFeatureCount.setForeground(Color.DARK_GRAY);
				checkFeatureCount.setBackground(Color.LIGHT_GRAY);
				checkFeatureCount.setBounds(10, 27, 138, 23);

				panelCounting.add(checkFeatureCount);
				checkHtSeq.setFont(GLSeq2_Main_Application.TEXT_FONT);
				checkHtSeq.setForeground(Color.DARK_GRAY);
				checkHtSeq.setBackground(Color.LIGHT_GRAY);
				checkHtSeq.setBounds(10, 48, 138, 23);

				panelCounting.add(checkHtSeq);
				checkRsem.setFont(GLSeq2_Main_Application.TEXT_FONT);
				checkRsem.setForeground(Color.DARK_GRAY);
				checkRsem.setBackground(Color.LIGHT_GRAY);
				checkRsem.setBounds(10, 69, 138, 23);

				panelCounting.add(checkRsem);
			}
			{
				txtchPipelineOptions.setBounds(10, 11, 284, 42);
				panel_1.add(txtchPipelineOptions);
				txtchPipelineOptions.setForeground(Color.DARK_GRAY);
				nimbusFix(Color.LIGHT_GRAY, txtchPipelineOptions);
				txtchPipelineOptions.setText("Pipeline Options");
				txtchPipelineOptions.setFont(GLSeq2_Main_Application.HEADER_FONT);
				txtchPipelineOptions.setEditable(false);
			}
			{
				panelNumberOfCores.setForeground(Color.DARK_GRAY);
				panelNumberOfCores.setBackground(Color.LIGHT_GRAY);
				panelNumberOfCores.setBounds(10, 277, 443, 40);
				panelManualHolder.add(panelNumberOfCores);
				panelNumberOfCores.setLayout(null);
				spinNumberOfCores.setModel(new SpinnerNumberModel(new Integer(4), new Integer(0), null, new Integer(1)));
				spinNumberOfCores.setBounds(380, 0, 53, 40);

				panelNumberOfCores.add(spinNumberOfCores);
				{
					txtcNumberOfCores.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtcNumberOfCores);
					txtcNumberOfCores.setEditable(false);
					txtcNumberOfCores.setBounds(0, 0, 370, 40);
					panelNumberOfCores.add(txtcNumberOfCores);
					txtcNumberOfCores.setText("Number of Cores to Use");
				}
			}
			{
				panelParallelExpressionComp.setForeground(Color.DARK_GRAY);
				panelParallelExpressionComp.setBackground(Color.LIGHT_GRAY);
				panelParallelExpressionComp.setLayout(null);
				panelParallelExpressionComp.setBounds(550, 277, 443, 40);
				panelManualHolder.add(panelParallelExpressionComp);
				{
					spinParallelExpression.setModel(new SpinnerNumberModel(new Integer(1), new Integer(0), null, new Integer(1)));
					spinParallelExpression.setBounds(380, 0, 53, 40);
					panelParallelExpressionComp.add(spinParallelExpression);
				}
				{
					txtcParallelComputationStreams
							.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtcParallelComputationStreams);
					txtcParallelComputationStreams.setEditable(false);
					txtcParallelComputationStreams
							.setText("Parallel Computation Streams for Expression Computation");
					txtcParallelComputationStreams.setBounds(0, 0, 370, 40);
					panelParallelExpressionComp.add(txtcParallelComputationStreams);
				}
			}
			{
				panelParallelDataComp.setForeground(Color.DARK_GRAY);
				panelParallelDataComp.setBackground(Color.LIGHT_GRAY);
				panelParallelDataComp.setLayout(null);
				panelParallelDataComp.setBounds(10, 328, 443, 40);
				panelManualHolder.add(panelParallelDataComp);
				{
					spinParallelDataPrep.setModel(new SpinnerNumberModel(new Integer(1), new Integer(0), null, new Integer(1)));
					spinParallelDataPrep.setBounds(380, 0, 53, 40);
					panelParallelDataComp.add(spinParallelDataPrep);
				}
				{
					txtcParallelComputationStreams_1
							.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY,
							txtcParallelComputationStreams_1);
					txtcParallelComputationStreams_1.setEditable(false);
					txtcParallelComputationStreams_1
							.setText("Parallel Computation Streams for Data Preparation");
					txtcParallelComputationStreams_1.setBounds(0, 0, 370, 40);
					panelParallelDataComp.add(txtcParallelComputationStreams_1);
				}
			}
			{
				panelMaxFragLen.setForeground(Color.DARK_GRAY);
				panelMaxFragLen.setBackground(Color.LIGHT_GRAY);
				panelMaxFragLen.setLayout(null);
				panelMaxFragLen.setBounds(550, 328, 443, 40);
				panelManualHolder.add(panelMaxFragLen);
				{
					spinMaxFragLen.setModel(new SpinnerNumberModel(new Integer(1000), new Integer(0), null, new Integer(100)));
					spinMaxFragLen.setBounds(380, 0, 53, 40);
					panelMaxFragLen.add(spinMaxFragLen);
				}
				{
					txtcMaximalFragmentLength.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtcMaximalFragmentLength);
					txtcMaximalFragmentLength.setEditable(false);
					txtcMaximalFragmentLength
							.setText("Maximal Fragment Length");
					txtcMaximalFragmentLength.setBounds(0, 0, 370, 40);
					panelMaxFragLen.add(txtcMaximalFragmentLength);
				}
			}
			{
				panelMaxBufferConfInts.setForeground(Color.DARK_GRAY);
				panelMaxBufferConfInts.setBackground(Color.LIGHT_GRAY);
				panelMaxBufferConfInts.setLayout(null);
				panelMaxBufferConfInts.setBounds(10, 384, 443, 40);
				panelManualHolder.add(panelMaxBufferConfInts);
				{
					spinMaxBuffer.setModel(new SpinnerNumberModel(new Integer(4096), new Integer(0), null, new Integer(1024)));
					spinMaxBuffer.setBounds(379, 0, 54, 40);
					panelMaxBufferConfInts.add(spinMaxBuffer);
				}
				{
					txtcMaximumAuxiliaryBuffer.setForeground(Color.DARK_GRAY);
					nimbusFix(Color.LIGHT_GRAY, txtcMaximumAuxiliaryBuffer);
					txtcMaximumAuxiliaryBuffer.setEditable(false);
					txtcMaximumAuxiliaryBuffer
							.setText("Maximum Auxiliary Buffer for Computing Credibility Intervals");
					txtcMaximumAuxiliaryBuffer.setBounds(0, 0, 370, 40);
					panelMaxBufferConfInts.add(txtcMaximumAuxiliaryBuffer);
				}
			}
			{
				btnExtractCoverage.setFont(GLSeq2_Main_Application.TEXT_FONT);
				btnExtractCoverage.setBounds(10, 229, 983, 33);
				btnExtractCoverage.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						if (btnExtractCoverage
								.getText()
								.equals(ButtonEnums.OptionButton.EXTRACT.value)) {
							btnExtractCoverage
									.setText(ButtonEnums.OptionButton.NO_EXTRACT.value);
						} else {
							btnExtractCoverage
									.setText(ButtonEnums.OptionButton.EXTRACT.value);
						}
					}
				});
				panelManualHolder.add(btnExtractCoverage);
			}
			{
				btnComputeConfIntervals.setFont(GLSeq2_Main_Application.TEXT_FONT);
				btnComputeConfIntervals.setBounds(10, 177, 443, 40);
				btnComputeConfIntervals.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						if (btnComputeConfIntervals.getText().equals(
								ButtonEnums.OptionButton.COMPUTE.value)) {
							btnComputeConfIntervals
									.setText(ButtonEnums.OptionButton.NO_COMPUTE.value);
						} else {
							btnComputeConfIntervals
									.setText(ButtonEnums.OptionButton.COMPUTE.value);
						}
					}
				});
				panelManualHolder.add(btnComputeConfIntervals);
			}
			{
				btnOutputGenomeBam.setFont(GLSeq2_Main_Application.TEXT_FONT);
				btnOutputGenomeBam.setBounds(550, 178, 443, 40);
				btnOutputGenomeBam.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						if (btnOutputGenomeBam.getText().equals(
								ButtonEnums.OptionButton.OUTPUT.value)) {
							btnOutputGenomeBam
									.setText(ButtonEnums.OptionButton.NO_OUTPUT.value);
						} else {
							btnOutputGenomeBam.setText(ButtonEnums.OptionButton.OUTPUT.value);
						}
					}
				});
				panelManualHolder.add(btnOutputGenomeBam);
			}
		}
		{
			buttonPane.setBackground(Color.GRAY);
			buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			{
				okButton.setActionCommand("OK");
				okButton.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent arg0) {
						GLSeq2_Main_Application.att.setaAlgor(String.valueOf(comboAlignmentAlgo.getSelectedItem()));
						GLSeq2_Main_Application.att.setRGenome(txtReferenceGenome.getText());
						GLSeq2_Main_Application.att.setFASTAname(txtReferenceFasta.getText());
						GLSeq2_Main_Application.att.setGFFname(txtReferenceFeatures.getText());
						GLSeq2_Main_Application.att.setIdAttr(txtReferenceID.getText());
						GLSeq2_Main_Application.att.setFeatureColumn(String.valueOf(spinColumnGtf.getValue()));
						if (String.valueOf(comboAlignmentAlgo.getSelectedItem()).equals("Cushaw-GPU")){
							GLSeq2_Main_Application.att.setGpuAccel("TRUE");
						} else{
							GLSeq2_Main_Application.att.setGpuAccel("FALSE");
						}
						GLSeq2_Main_Application.att.setCores(String.valueOf(spinNumberOfCores.getValue()));
						if (checkFeatureCount.isSelected()){
							GLSeq2_Main_Application.att.setFeatureCounts("FeatureCounts");
						} else{
							GLSeq2_Main_Application.att.setFeatureCounts("");
						}
						if (checkRsem.isSelected()){
							GLSeq2_Main_Application.att.setRSEM("RSEM");
						} else{
							GLSeq2_Main_Application.att.setRSEM("");
						}
						if (checkHtSeq.isSelected()){
							GLSeq2_Main_Application.att.setHTSeq("HTSeq");
						} else{
							GLSeq2_Main_Application.att.setHTSeq("");
						}
						GLSeq2_Main_Application.run.setProtocolId(txtProtocolID.getText());
						GLSeq2_Main_Application.att.setStreams(String.valueOf(spinParallelExpression.getValue()));
						GLSeq2_Main_Application.att.setStreamsDataPrep(String.valueOf(spinParallelDataPrep.getValue()));
						GLSeq2_Main_Application.att.setFragMaxLength(String.valueOf(spinMaxFragLen.getValue()));
						GLSeq2_Main_Application.att.setCiMem(String.valueOf(spinMaxBuffer.getValue()));
						
						if (btnExtractCoverage.getText().equals(ButtonEnums.OptionButton.EXTRACT.value)){
							GLSeq2_Main_Application.att.setStrandExtract("TRUE");
						} else{
							GLSeq2_Main_Application.att.setStrandExtract("FALSE");
						}
						
						if (btnComputeConfIntervals.getText().equals(ButtonEnums.OptionButton.COMPUTE.value)){
							GLSeq2_Main_Application.att.setCompConf("TRUE");
						} else{
							GLSeq2_Main_Application.att.setCompConf("FALSE");
						}
						
						if (btnOutputGenomeBam.getText().equals(ButtonEnums.OptionButton.OUTPUT.value)){
							GLSeq2_Main_Application.att.setGenoBam("TRUE");
						} else{
							GLSeq2_Main_Application.att.setGenoBam("FALSE");
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
	void nimbusFix(Color background, JTextPane pane) {
		UIDefaults defaults = new UIDefaults();
		defaults.put("TextPane[Enabled].backgroundPainter", background);
		pane.putClientProperty("Nimbus.Overrides", defaults);
		pane.putClientProperty("Nimbus.Overrides.InheritDefaults", true);
		pane.setBackground(background);
	}
}
