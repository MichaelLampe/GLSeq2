package application;

import java.util.ArrayList;
import java.util.List;

public class Run {
  // Static keep count stuff
  private static int runs_count = 0;
  private RunTableEntry active_run;
  // Table Data
  public final String run_name;
  private final boolean dataprep;
  private final boolean alignment;
  private final boolean counting;
  private final boolean collect;

  public Run(boolean dataprep, boolean alignment, boolean counting, boolean collect, String run_name) {
    this.run_name = run_name;
    this.dataprep = dataprep;
    this.alignment = alignment;
    this.counting = counting;
    this.collect = collect;
    runs_count = runs_count + 1;
    active_run = new RunTableEntry(this.run_name);
  }

  public static int getActiveRunCount() {
    return runs_count;
  }

  public RunTableEntry getActiveRun() {
    return active_run;
  }

  public List<String> constructArgs() {
    List<String> args = new ArrayList<String>();
    args.add("Rscript");
    args.add("GLSeq.top.R");
    args.add("Placeholder");
    if (dataprep) {
      args.add("dataprep");
    } else {
      args.add("nodataprep");
    }
    if (alignment) {
      args.add("alignment");
    } else {
      args.add("noalignment");
    }
    if (counting) {

      args.add("counting");
    } else {
      args.add("nocounting");
    }
    if (collect) {
      args.add("collect");
    } else {
      args.add("nocollect");
    }
    args.add(run_name);
    args.add("0");
    return args;
  }
}
