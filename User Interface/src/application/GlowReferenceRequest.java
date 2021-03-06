package application;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.jcabi.ssh.Shell;

import javafx.application.Platform;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;

public class GlowReferenceRequest extends GlowRequest {

	private String reference_id;
	/*
	 * Names of Attribute hashes to assign values.
	 */
	private static final String REFERENCE_FASTA_NAME = "referenceFasta";
	private static final String REFERENCE_GFF_NAME = "referenceGff";
	private static final String REFERENCE_GENOME_NAME = "referenceGenome";

	public GlowReferenceRequest(GlowRequest connection, String reference_id) {
		super(connection.getUsername(), connection.getPassword());
		this.reference_id = reference_id;
		this.connected = connection.getConnected();
	}

	private List<String> requestReferenceFiles(String reference_id) {
		/*
		 * Simple curl command to ask for stuff.
		 * 
		 * TODO: Add in SSH library again to allow for this to be run over an
		 * SSH.
		 */
		@SuppressWarnings("serial")
		ArrayList<String> command = new ArrayList<String>() {
			{
				add("curl");
				add("--cookie");
				add(cookie);
				add(glow_server + "reference_genomes/" + reference_id);
			}
		};
		return command;
	}

	private List<String> requestReferenceName(String reference_id) {
		@SuppressWarnings("serial")
		/*
		 * Simple curl command to ask for stuff.
		 * 
		 * TODO: Add in SSH library again to allow for this to be run over an
		 * SSH.
		 */
		ArrayList<String> command = new ArrayList<String>() {
			{
				add("curl");
				add("--cookie");
				add(cookie);
				add(glow_server + "reference_genomes/" + reference_id + ".json");
			}
		};
		return command;
	}

	@Override
	protected Object call() throws Exception {
		/*
		 * Cookie request takes a while, so give the alert prior to give a
		 * notice to the user. Needs to run on platform thread so we call it as
		 * a lambda.
		 */
		Platform.runLater(() -> {
			Alert alertLooking = new Alert(AlertType.INFORMATION);
			alertLooking.setTitle("Reference Files");
			alertLooking.setHeaderText("Requesting Reference Files from GLOW Server");
			alertLooking.setContentText("Another alert will occur when the files have been received.");
			alertLooking.show();
		});

		/*
		 * Check the cookie, if it is not good try to renew it
		 */
		if (renewCookieExpiration()) {

			/*
			 * We make two requests here, one for files and one for the name.
			 * The reference name is gotten from a .json request of the same
			 * page as the reference files request.
			 */

			BufferedReader fileR = null;
			BufferedReader nameR = null;

			/*
			 * We can parse input streams to only add what we want to our list.
			 * 
			 * fileR: This pretty much says that if a lien contains either of
			 * the two strings below, add them to our commands list.
			 * 
			 * nameR: The same is true for the this, if we find a line that
			 * contains name we add it to our command line.
			 */
			List<String> commands = new ArrayList<String>();
			// Linux version
			if (runningOnLinux) {
				Process processFiles = new ProcessBuilder(requestReferenceFiles(reference_id)).start();
				Process processName = new ProcessBuilder(requestReferenceName(reference_id)).start();

				/*
				 * Record the values returned via an input stream.
				 */
				fileR = new BufferedReader(new InputStreamReader(processFiles.getInputStream()));
				nameR = new BufferedReader(new InputStreamReader(processName.getInputStream()));

				/*
				 * Wait for streams to finish.
				 */
				processFiles.waitFor();
				processName.waitFor();

			} else {
				// Windows
				Shell ssh = establishSsh();
				try {
					System.out.println("Looking for files");
					fileR = new BufferedReader(new StringReader(new Shell.Plain(ssh)
							.exec(convertCommandToString((ArrayList<String>) requestReferenceFiles(reference_id)))));
					nameR = new BufferedReader(new StringReader(new Shell.Plain(ssh)
							.exec(convertCommandToString((ArrayList<String>) requestReferenceName(reference_id)))));
				} catch (Exception e) {
					e.printStackTrace();
				}
			}

			commands.addAll(fileR.lines()
					.filter(line -> (line.contains("<td>GFF</td>")) || (line.contains("<td>Sequence:FastA</td>")))
					.collect(Collectors.toList()));

			commands.addAll(nameR.lines().map(line -> line.split(","))
					// Flat map makes it happy that we split it again.
					.flatMap(Arrays::stream).filter(line -> (line.contains("name"))).collect(Collectors.toList()));
			/*
			 * Update our Attributes based on the values we just had returned
			 * from GLOW. If this doesn't work, they may have changed the format
			 * of the response.
			 */
			Platform.runLater(() -> {
				Alert alertReceived = new Alert(AlertType.INFORMATION);
				alertReceived.setTitle("Received Reference Files");
				alertReceived.setHeaderText("Files have been received from GLOW.");
				alertReceived.setContentText("Fields populated with all data retrieved from GLOW for that ID");
				alertReceived.show();
			});
			for (String line : commands) {
				if (line.contains("name")) {
					/*
					 * Parsing
					 */
					String[] temp = line.split("\"");
					String referenceGenomeFileName = temp[3];

					// Updating
					AttributeFactorySingleton.getInstance().setAttributeValue(REFERENCE_GENOME_NAME,
							referenceGenomeFileName);
				} else if (line.contains("<td>GFF</td>")) {
					/*
					 * Parsing
					 */
					String[] temp = line.split("class=\"detail_link\">");
					String gffFileName = temp[1].split("</a>")[0];

					// Updating
					AttributeFactorySingleton.getInstance().setAttributeValue(REFERENCE_GFF_NAME, gffFileName);
				} else {
					/*
					 * Parsing
					 */
					String[] temp = line.split("class=\"detail_link\">");
					String fastaFileName = temp[1].split("</a>")[0];

					// Updating
					AttributeFactorySingleton.getInstance().setAttributeValue(REFERENCE_FASTA_NAME, fastaFileName);
				}
			}

			/*
			 * Only JavaFx thread can update JavaFx. We can call this thread to
			 * update via a Platform lambda that indicates the defaults of the
			 * UI should be updated.
			 */
			Platform.runLater(() -> {
				UpdateUserInterfaceSingleton.getInstance().updateDefaults();
			});
		} else {
			System.out.println("Did not run because no connection to glow was established.");
		}
		return null;
	}
}
