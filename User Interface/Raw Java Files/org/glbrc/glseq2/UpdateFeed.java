package org.glbrc.glseq2;

import java.util.Timer;
import java.util.TimerTask;

import javax.swing.JTextArea;

public class UpdateFeed extends JTextArea {
  private static final long serialVersionUID = 1L;
  private static String currentText = "";
  private Timer time;

  /**
   * Initializes the update feed class.
   */
  public UpdateFeed() {
    // Create at current time
    setWrapStyleWord(true);
    setLineWrap(true);
    time = new Timer();
    currentText = getText();
  }

  /**
   * Updates the text field.
   * 
   * @return true
   */
  public boolean update(String update) {
    // New lines create white space and improves readability
    updateText("\n" + update + "\n");
    //
    // Solution to the update problems. Found on stack overflow.
    // http://stackoverflow.com/a/11166024
    //
    time.schedule(new TimerTask() {
      @Override
      public void run() {
        try {
          setText(currentText);
          moveView();
        } catch (NullPointerException e) {
          // About 1 in 10 runs a null pointer will pop up here.
          // Decided to just add a print statement for debugging so the user
          // Won't see a huge stack trace, as the error doesn't seem to affect
          // function.
          System.out.println("Possible timer problem. UpdateFeed.java");
        }
      }
    }, 100);
    // Keep the scroll bar at the bottom
    return true;
  }

  /*
   * Moves the view to the bottom of the viewport
   */
  private void moveView() {
    setCaretPosition(this.getDocument().getLength());
  }

  /**
   * Updates the text field and removes extra characters. If there are too many.
   * 
   * @param the
   *          String to add
   */
  private static void updateText(String update) {
    currentText += update;
    while (currentText.length() > 5000) {
      currentText = currentText.substring(500, currentText.length());
    }
  }
}
