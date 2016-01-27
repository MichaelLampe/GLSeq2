package application;

import java.util.ArrayList;
import java.util.List;

import javafx.concurrent.Task;

public class RunInstantiator {
	// Static keep count stuff
	private static int runs_count = 0;
	// Table Data
	public final String run_name;
	private final boolean dataprep;
	private final boolean alignment;
	private final boolean counting;
	private final boolean collect;
	private List<String> args = new ArrayList<String>();

	public RunInstantiator(boolean dataprep, boolean alignment, boolean counting,
			boolean collect, String run_name) {
		this.run_name = run_name;
		this.dataprep = dataprep;
		this.alignment = alignment;
		this.counting = counting;
		this.collect = collect;
		runs_count = runs_count + 1;
	}

	public static int getActiveRunCount() {
		return runs_count;
	}

	public void constructArgs(boolean Condor, String attribute_file,
			String scriptsDirectory) {
		// Condor Wrapper to run via Condor
		if (!scriptsDirectory.endsWith("/")) {
			scriptsDirectory = scriptsDirectory + "/";
		}
		if (Condor) {
			args.add("python");
			args.add(scriptsDirectory + "PyGLSeqWrapper.py");
			// RScript
		} else {
			args.add("Rscript");
			args.add(scriptsDirectory + "GLSeq.top.R");
		}
		// Where update was, might add to here later
		args.add("Placeholder");

		// Dataprep
		if (dataprep)
			args.add("dataprep");
		else
			args.add("nodataprep");

		// Alignment
		if (alignment)
			args.add("alignment");
		else
			args.add("noalignment");

		// Counting
		if (counting)
			args.add("counting");
		else
			args.add("nocounting");

		// Collect
		if (collect)
			args.add("collect");
		else
			args.add("nocollect");

		// This means it is a directory, so we don't assign a name
		if (!attribute_file.endsWith("/")) {
			args.add(this.run_name);
		}

		// Protocol ID
		args.add("0");

		// Attribute File
		args.add(attribute_file);
	}

	public void start() {
		if (args != null) {
			/*
			 * Might add some logic in the future that starts an SSH on
			 * windows/mac
			 */
			Task<Object> task = new GlSeqRunWorker(args);
			// Run on a different thread so it doesn't lock up the UI.
			Thread thread = new Thread(task);
			thread.setDaemon(true);
			thread.start();
		} else {
			System.out.println("Please select something for the script to do.");
		}
	}
}
