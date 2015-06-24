package org.glbrc.glseq2;

import java.awt.Component;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JTabbedPane;

public class BatchTab extends JTabbedPane {

  private static final long serialVersionUID = 1L;

  public BatchTab(int right) {
    super(right);
  }

  // These two methods make it so only QueuedRun classes can be added
  // Nice wrappers that solves some problems that can come up
  // if we tried to constantly cast.
  public void addQueue(QueuedRun run) {
    add(run);
  }

  public void removeQueue(QueuedRun run) {
    remove(run);
  }

  /**
   * Forces the return of only QueuedRuns.
   * 
   * @return Returns all QueuedRuns attached to this pane.
   */
  public List<QueuedRun> getQueues() {
    List<QueuedRun> queuedRuns = new ArrayList<QueuedRun>();
    for (Component comp : getComponents()) {
      if (comp instanceof QueuedRun) {
        QueuedRun run = (QueuedRun) comp;
        queuedRuns.add(run);
      }
    }
    return queuedRuns;
  }
}
