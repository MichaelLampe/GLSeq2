package org.glbrc.glseq2;

import java.awt.Color;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;
import javax.swing.ScrollPaneConstants;

public class QueuedRun extends JPanel {
  public static int count = 0;

  private final Color onColor = new Color(0, 128, 0);
  private final Color offColor = new Color(178, 34, 34);
  private static final long serialVersionUID = 1L;
  private final JButton btnUpdate = new JButton("Update");
  private final JButton btnDataprep = new JButton("DataPrep");
  private final JButton btnAlignment = new JButton("Alignment");
  private final JButton btnCounting = new JButton("Counting");
  private final JButton btnCollecting = new JButton("Collecting");
  private final JTextPane txtAttributePath = new JTextPane();
  private final JTextPane txtRunName = new JTextPane();
  private final JTextPane txtAlignment = new JTextPane();
  private final JTextPane txtCounting = new JTextPane();
  private final JTextPane txtPaired = new JTextPane();
  private final JTextPane txtZipped = new JTextPane();

  private final Attributes panelAttributes;
  private final Run panelRun;
  private final JScrollPane scrollPane = new JScrollPane();
  private final JScrollPane scrollPane1 = new JScrollPane();
  private final JButton deleteRun = new JButton("X");

  /**
   * Create the panel.
   */
  public QueuedRun(Run panelRun, Attributes panelAttributes) {
    this.panelRun = panelRun;
    this.panelAttributes = panelAttributes;
    count++;
    initGui();
  }

  private void initGui() {
    setLayout(null);
    btnUpdate.setForeground(Color.WHITE);
    btnUpdate.setBounds(170, 36, 142, 23);
    if (panelRun.getUpdateFromDatabase().equals(ButtonEnums.Attribute.UPDATE.value)) {
      btnUpdate.setBackground(onColor);
      btnUpdate.setText("Updating");
    } else {
      btnUpdate.setBackground(offColor);
      btnUpdate.setText("Not Updating");
    }
    add(btnUpdate);
    btnDataprep.setForeground(Color.WHITE);
    btnDataprep.setBounds(170, 66, 142, 23);
    if (panelRun.getProcessedData().equals(ButtonEnums.Attribute.PREPROCESSING.value)) {
      btnDataprep.setText("Preparing Data");
      btnDataprep.setBackground(onColor);
    } else {
      btnDataprep.setBackground(offColor);
      btnDataprep.setText("Not Preparing Data");
    }

    add(btnDataprep);
    btnAlignment.setForeground(Color.WHITE);
    btnAlignment.setBounds(170, 96, 142, 23);
    if (panelRun.getAlignment().equals(ButtonEnums.Attribute.ALIGNMENT.value)) {
      btnAlignment.setText("Aligning");
      btnAlignment.setBackground(onColor);
    } else {
      btnAlignment.setBackground(offColor);
      btnAlignment.setText("Not Aligning");
    }

    add(btnAlignment);
    btnCounting.setForeground(Color.WHITE);
    btnCounting.setBounds(170, 126, 142, 23);
    if (panelRun.getCounting().equals(ButtonEnums.Attribute.COUNT.value)) {
      btnCounting.setText("Counting");
      btnCounting.setBackground(onColor);
    } else {
      btnCounting.setBackground(offColor);
      btnCounting.setText("Not Counting");
    }

    add(btnCounting);
    btnCollecting.setForeground(Color.WHITE);
    btnCollecting.setBounds(170, 156, 142, 23);
    if (panelRun.getCollectResults().equals(ButtonEnums.Attribute.COLLECT.value)) {
      btnCollecting.setText("Collecting");
      btnCollecting.setBackground(onColor);
    } else {
      btnCollecting.setBackground(offColor);
      btnCollecting.setText("Not Collecting");
    }

    add(btnCollecting);
    txtRunName.setEditable(false);
    txtRunName.setText(panelRun.getRunId());
    txtRunName.setBounds(10, 10, 250, 20);

    add(txtRunName);
    txtAlignment.setEditable(false);
    if (panelAttributes.getaAlgor().equals("Cushaw")) {
      if (panelAttributes.getGpuAccel().equals("TRUE")) {
        txtAlignment.setText("Cushaw-GPU");
      } else {
        txtAlignment.setText(panelAttributes.getaAlgor());
      }
    } else {
      txtAlignment.setText(panelAttributes.getaAlgor());
    }

    txtAlignment.setBounds(10, 40, 150, 20);

    add(txtAlignment);
    if (panelAttributes.getPairedEnd().equals("TRUE")) {
      txtPaired.setText(ButtonEnums.OptionButton.PAIRED.value);
    } else {
      txtPaired.setText(ButtonEnums.OptionButton.SINGLE.value);
    }
    txtPaired.setEditable(false);

    txtPaired.setBounds(10, 130, 150, 20);

    add(txtPaired);
    if (panelAttributes.getUnzipped().equals("TRUE")) {
      txtZipped.setText(ButtonEnums.OptionButton.UNZIPPED.value);
    } else {
      txtZipped.setText(ButtonEnums.OptionButton.ZIPPED.value);
    }
    txtZipped.setEditable(false);

    txtZipped.setBounds(10, 160, 150, 20);

    add(txtZipped);
    scrollPane.setBounds(10, 191, 302, 44);

    add(scrollPane);
    txtAttributePath.setEditable(false);
    scrollPane.setViewportView(txtAttributePath);
    txtAttributePath.setText(panelRun.getAttributeFilePath());
    scrollPane1.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
    scrollPane1.setBounds(10, 69, 150, 50);

    add(scrollPane1);
    scrollPane1.setViewportView(txtCounting);
    txtCounting.setEditable(false);
    txtCounting.setText(countingString());
    deleteRun.setForeground(Color.RED);
    deleteRun.setFont(new Font("Arial Black", Font.BOLD, 11));
    deleteRun.setBounds(269, 10, 43, 20);
    deleteRun.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent btn) {
        removeThis();
      }
    });
    add(deleteRun);
  }

  /**
   * Returns a string that a human would like to read of the countable variables
   * being used.
   * 
   * @return countString, in format for human to read
   */
  private String countingString() {
    String countString = "";
    if (panelAttributes.getRsem().equals("RSEM")) {
      countString += panelAttributes.getRsem();
    }
    if (panelAttributes.getFeatureCounts().equals("FeatureCounts")) {
      if (!countString.equals("")) {
        countString += "\n";
      }
      countString += panelAttributes.getFeatureCounts();
    }
    if (panelAttributes.getHtseq().equals("HTSeq")) {
      if (!countString.equals("")) {
        countString += "\n";
      }
      countString += panelAttributes.getHtseq();
    }
    if (panelAttributes.getCufflinks().equals("Cufflinks")) {
      if (!countString.equals("")) {
        countString += "\n";
      }
      countString += panelAttributes.getCufflinks();
    }

    return countString;
  }

  private void removeThis() {
    Application.tabsRun.remove(this);
    count--;
  }

  public Run getSelectedRun() {
    return panelRun;
  }

  public Attributes getSelectedAttributes() {
    return panelAttributes;
  }
}
