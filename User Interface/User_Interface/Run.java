package application;

import java.util.ArrayList;
import java.util.List;

public class Run {
	// Static keep count stuff
	private static int runs_count = 0;
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
	}

	public static int getActiveRunCount() {
		return runs_count;
	}

	public List<String> constructArgs(boolean Condor, String attribute_file) {
		List<String> args = new ArrayList<String>();
		// Condor Wrapper to run via Condor
		if (Condor) {
			args.add("python");
			args.add("PyGLSeqWrapper.py");
			// RScript
		} else {
			args.add("Rscript");
			args.add("GLSeq.top.R");
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
		if (!attribute_file.endsWith("/"))
			args.add(run_name);

		// Protocol ID
		args.add("0");

		// Attribute File
		args.add(attribute_file);

		return args;
	}
}
