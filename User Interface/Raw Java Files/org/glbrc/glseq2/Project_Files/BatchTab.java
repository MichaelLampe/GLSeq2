package org.glbrc.glseq2.Project_Files;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JTabbedPane;

public class BatchTab extends JTabbedPane {

  private static final long serialVersionUID = 1L;

  private List<QueuedRun> queues;

  public BatchTab(int right) {
    super(right);
    queues = new ArrayList<QueuedRun>();

  }

  // These two methods make it so only QueuedRun classes can be added
  // Nice wrappers that solves some problems that can come up
  // if we tried to constantly cast.
  /**
   * Adds the run to the queue.
   * 
   * @param run
   *          - The run added to the queue
   */
  public void addQueue(QueuedRun run) {
    add(run);
    queues.add(run);
  }

  public void removeQueue(QueuedRun run) {
    remove(run);
    queues.remove(run);
  }

  /**
   * Forces the return of only QueuedRuns.
   * 
   * @return Returns all QueuedRuns attached to this pane.
   */
  public List<QueuedRun> getQueues() {
    return getQueue_NoClientCode();
  }

  /**
   * This is the structure utilized by the main JTabbedPane.
   */
  private List<QueuedRun> getQueue_NoClientCode() {
    synchronized (getTreeLock()) {
      return queues;
    }
  }
/**
 * Removes all the queues from this tab.
 */
  public void removeQueues() {
    // Removes all the used tabs.
    int remove = queues.size() - 1;
    for (int i = remove; i >= 0; i--) {
      Application.tabsRun.removeQueue(queues.get(i));
      QueuedRun.count--;
    }
  }
}
