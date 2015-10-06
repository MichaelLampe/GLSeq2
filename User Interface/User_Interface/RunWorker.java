package application;

import java.io.File;
import java.util.List;

import javafx.concurrent.Task;

public class RunWorker extends Task<Object> {

	private List<String> args;

	public RunWorker(List<String> args) {
		this.args = args;
	}

	@Override
	protected Object call() throws Exception {
		// Script generated based on arguments from the RunOptions class.
		// It calls the R script with the correct user parameters
		// Build process w/ args again
		String scriptDirectory = Attributes.getInstance().attributesCollection.get("base.dir").getValue();
		Process process = null;
		try {
			process = new ProcessBuilder(args).directory(new File(scriptDirectory)).inheritIO().start();
			// Wait until it is done or it crashes.
			process.waitFor();
		} catch (Exception e) {
			System.out.println("Error in the process.");
			process.destroy();
		}
		return null;
	}
}
