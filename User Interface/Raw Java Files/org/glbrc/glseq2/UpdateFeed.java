package org.glbrc.glseq2;

import java.util.Date;
import java.util.Timer;
import java.util.TimerTask;

import javax.swing.JTextArea;

public class UpdateFeed extends JTextArea {
  private static final long serialVersionUID = 1L;

  /**
   * Create the panel.
   */
  private static Date currentTime;
  private static Date lastUpdate;
  private static String currentText;
  private static boolean futureRun = false;
  private Timer time;

  /**
   * Initializes the update feed class.
   */
  public UpdateFeed() {
    // Create at current time
    setWrapStyleWord(true);
    setLineWrap(true);
    time = new Timer();
    lastUpdate = new Date();
    currentText = getText();
  }

  /**
   * Updates the text field.
   * 
   * @return true
   */
  public boolean update(String update) {
    if (futureRun) {
      // Yes I did mean to add to new lines, this creates a nice white space
      // around individual messages.
      updateText("\n" + update + "\n");
      futureRun = false;
    }
    if (checkTime()) {
      updateText("\n" + update + "\n");
      setText(currentText);
    } else {
      // Solution to the update problems. Found on stack overflow.
      // http://stackoverflow.com/a/11166024
      //
      futureRun = true;
      time.schedule(new TimerTask() {
        @Override
        public void run() {
          setText(currentText);
        }
      }, 100);
    }
    // Keep the scrollbar at the bottom
    setCaretPosition(this.getDocument().getLength());
    return true;
  }

  /**
   * Makes the text field wait 2 seconds to update.
   * 
   * @return True if waiting is no longer needed
   */
  private static boolean checkTime() {
    currentTime = new Date();
    if ((currentTime.getTime() - lastUpdate.getTime()) > 100) {
      lastUpdate = new Date();
      return true;
    } else {
      return false;
    }
  }

  private static void updateText(String update) {
    currentText += update;
    while (currentText.length() > 5000) {
      currentText = currentText.substring(500, currentText.length());
    }
  }
}
