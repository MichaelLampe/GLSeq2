package org.glbrc.glseq2;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.JTextPane;
import javax.swing.SpinnerNumberModel;
import javax.swing.UIDefaults;
import javax.swing.border.EmptyBorder;

public class Processing extends JDialog {

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
  private final JTextPane txtcMinimumTrimRead = new JTextPane();
  private final JTextPane txtchTrimmingAndProcessing = new JTextPane();
  private final JPanel buttonPane = new JPanel();
  private final JButton okButton = new JButton("Apply and Close");
  private final JButton cancelButton = new JButton("Cancel");

  /**
   * Launch the application.
   * 
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
   * 
   */
  public Processing() {
    initGui();
    // Assign some variables
    txtArtificialName.setText(Application.att.getArtificialFasta());
    spinTrimHead.setValue(Integer.valueOf(Application.att.getTrimhead()));
    spinMinTrimRead.setBounds(473, 0, 98, 40);
    panelTrimHead.add(spinMinTrimRead);
    spinMinTrimRead.setValue(Integer.valueOf(Application.att.getMinTrim()));
    if (Application.att.getReadTrim().equals("FALSE")) {
      btnTrimReads.setText(ButtonEnums.OptionButton.NO_TRIMMING.value);
    }

  }

  private void initGui() {
    setBackground(Color.LIGHT_GRAY);
    setResizable(false);
    setBounds(100, 100, 818, 326);
    getContentPane().setLayout(new BorderLayout());
    contentPanel.setBackground(Color.LIGHT_GRAY);
    contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    getContentPane().add(contentPanel, BorderLayout.CENTER);
    contentPanel.setLayout(null);
    {
      panelArtificialSequence.setForeground(Color.DARK_GRAY);
      panelArtificialSequence.setBackground(Color.LIGHT_GRAY);
      panelArtificialSequence.setLayout(null);
      panelArtificialSequence.setBounds(10, 131, 792, 59);
      contentPanel.add(panelArtificialSequence);
      {
        txtArtificialName.setBounds(298, 11, 484, 37);
        panelArtificialSequence.add(txtArtificialName);
      }
      {
        txtcNameOfFasta.setForeground(Color.DARK_GRAY);
        nimbusFix(Color.LIGHT_GRAY, txtcNameOfFasta);
        txtcNameOfFasta.setText("Name of FASTA File with Artificial Sequences");
        txtcNameOfFasta.setFont(Application.TEXT_FONT);
        txtcNameOfFasta.setEditable(false);
        txtcNameOfFasta.setBounds(10, 11, 278, 37);
        panelArtificialSequence.add(txtcNameOfFasta);
      }
    }
    {
      panelTrimHead.setForeground(Color.DARK_GRAY);
      panelTrimHead.setBackground(Color.LIGHT_GRAY);
      panelTrimHead.setLayout(null);
      panelTrimHead.setBounds(10, 201, 792, 40);
      contentPanel.add(panelTrimHead);
      {
        spinTrimHead.setModel(new SpinnerNumberModel(new Integer(12), new Integer(0), null,
            new Integer(1)));
        spinTrimHead.setBounds(175, 0, 98, 40);
        panelTrimHead.add(spinTrimHead);
      }
      {
        txtcTrimHead.setForeground(Color.DARK_GRAY);
        nimbusFix(Color.LIGHT_GRAY, txtcTrimHead);
        txtcTrimHead.setText("Trim Head");
        txtcTrimHead.setEditable(false);
        txtcTrimHead.setBounds(0, 0, 165, 40);
        panelTrimHead.add(txtcTrimHead);
      }
    }
    {
      {
        spinMinTrimRead.setModel(new SpinnerNumberModel(new Integer(36), new Integer(0), null,
            new Integer(1)));
      }
      {
        txtcMinimumTrimRead.setBounds(298, 0, 165, 40);
        panelTrimHead.add(txtcMinimumTrimRead);
        txtcMinimumTrimRead.setForeground(Color.DARK_GRAY);
        nimbusFix(Color.LIGHT_GRAY, txtcMinimumTrimRead);
        txtcMinimumTrimRead.setText("Minimum Trim Read");
        txtcMinimumTrimRead.setEditable(false);
      }
    }
    {
      btnTrimReads.setForeground(Color.DARK_GRAY);
      btnTrimReads.setFont(new Font("Arial", Font.PLAIN, 20));
      btnTrimReads.setBounds(10, 51, 792, 69);
      btnTrimReads.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent arg0) {
          if (btnTrimReads.getText().equals(ButtonEnums.OptionButton.TRIMMING.value)) {
            btnTrimReads.setText(ButtonEnums.OptionButton.NO_TRIMMING.value);
          } else {
            btnTrimReads.setText(ButtonEnums.OptionButton.TRIMMING.value);
          }
        }
      });
      contentPanel.add(btnTrimReads);
    }
    {
      nimbusFix(Color.LIGHT_GRAY, txtchTrimmingAndProcessing);
      txtchTrimmingAndProcessing.setForeground(Color.DARK_GRAY);
      txtchTrimmingAndProcessing.setFont(Application.HEADER_FONT);
      txtchTrimmingAndProcessing.setText("Trimming and Processing Options");
      txtchTrimmingAndProcessing.setEditable(false);
      txtchTrimmingAndProcessing.setBounds(10, 7, 792, 33);
      contentPanel.add(txtchTrimmingAndProcessing);
    }
    {
      buttonPane.setBackground(Color.GRAY);
      buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
      getContentPane().add(buttonPane, BorderLayout.SOUTH);
      {
        okButton.setActionCommand("OK");
        okButton.addActionListener(new ActionListener() {
          public void actionPerformed(ActionEvent arg0) {
            Application.att.setArtificialFasta(txtArtificialName.getText());
            Application.att.setMinTrim(String.valueOf(spinMinTrimRead.getValue()));
            Application.att.setTrimhead(String.valueOf(spinTrimHead.getValue()));
            if (btnTrimReads.getText().equals(ButtonEnums.OptionButton.NO_TRIMMING.value)) {
              Application.att.setReadTrim("FALSE");
            } else {
              Application.att.setReadTrim("TRUE");
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
