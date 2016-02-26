package application;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Locale;

import com.jcabi.ssh.SSHByPassword;
import com.jcabi.ssh.Shell;

import javafx.concurrent.Task;

/*
 * Glow request handles all the login stuff so that classes that extend it can handle the further details of their individual request implementations
 */
public class GlowRequest extends Task<Object> {
	protected final String glow_server = "https://glow-trunk.glbrc.org/";
	// Condor GLBRC "Scarcity-cm.glbrc.org" ip address
	protected final String sshServer = "144.92.98.39";

	protected String username;
	/*
	 * While a string is not the best way to store this, the password Will be
	 * logged in any string we make that calls a curl command anyway.
	 */
	protected String password = "";
	protected String cookie;
	protected boolean runningOnLinux;
	// If connected to GLOW by cookie and ping was true
	protected boolean connected = false;
	protected Shell localShell = null;

	public GlowRequest(String username, String password) {
		this.username = username;
		this.password = password;
		this.cookie = "cjar_" + username;
		this.runningOnLinux = true;

		String operatingSystem = System.getProperty("os.name", "generic").toLowerCase(Locale.ENGLISH);

		/*
		 * Not running on linux
		 */
		if (operatingSystem.indexOf("nux") < 0) {
			this.runningOnLinux = false;
		}
	}

	public Shell establishSsh() {
		// Return already established shell
		if (localShell != null) {
			return localShell;
		}

		try {
			Shell shell = new SSHByPassword(sshServer, 22, this.username, this.password);

			// Test connection
			new Shell.Plain(shell).exec("echo 'test'");

			localShell = shell;
			return shell;
		} catch (IllegalArgumentException | IOException e) {
			// Comes up if we try to log off after logging out.
			return null;
		}
	}

	String getUsername() {
		return username;
	}

	String getPassword() {
		return password;
	}

	boolean getConnected() {
		return connected;
	}

	String convertCommandToString(ArrayList<String> commandArray) {
		StringBuilder command = new StringBuilder("");
		for (int i = 0; i < commandArray.size(); i++) {
			if (i > 0) {
				command.append(" ");
			}
			command.append(commandArray.get(i));
		}
		return command.toString();
	}

	boolean requestCookie() {
		if (runningOnLinux) {
			ArrayList<String> command = new ArrayList<String>() {
				/**
				 * 
				 */
				private static final long serialVersionUID = 1L;

				{
					add("curl");
					add("--cookie-jar");
					add(cookie);
					add("--data");
					add("user[login]=" + username + "&user[password]=" + password + "&commit=Sign+in");
					add(glow_server + "users/sign_in");
				}
			};
			try {
				String response = "";
				Process process = new ProcessBuilder(command).directory(new File(System.getProperty("user.dir")))
						.start();
				BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(process.getInputStream()));
				String line = null;
				while ((line = bufferedReader.readLine()) != null) {
					response += line;
				}
				process.waitFor();
				/*
				 * The actual response if correct should be about 100, but this
				 * future proofs it a bit. If the password is wrong it returns
				 * the whole html page so even 200 should be pretty safe. The
				 * greater than protects us against incorrectly directed
				 * requests that don't do anything Nor get anything back making
				 * the program think they are correct
				 */
				return (response.length() < 200 && response.length() > 2);
			} catch (IOException | InterruptedException e) {
				/*
				 * Don't print stack trace here as it would reveal the user's
				 * password in that trace
				 */
			}
			return false;
		} else {
			Shell ssh = establishSsh();
			ArrayList<String> command = new ArrayList<String>() {

				/**
				 * 
				 */
				private static final long serialVersionUID = 1L;

				{
					add("curl");
					add("--cookie-jar");
					add(cookie);
					add("--data");
					add("'user[login]=" + username + "&user[password]=" + password + "&commit=Sign+in'");
					add(getLoginSite());
				}
			};
			try {
				// This keeps the request off the main draw thread and also
				// makes the request go a lot faster, win win.
				new Thread(new Shell.Plain(ssh).exec(convertCommandToString(command))).start();
			} catch (

			Exception e)

			{
				return false;
			}
			return true;

		}
	}

	String getLoginSite() {
		return glow_server + "users/sign_in";
	}

	boolean pingGlow() {
		@SuppressWarnings("serial")
		ArrayList<String> command = new ArrayList<String>() {
			{
				add("curl");
				add("--cookie");
				add(cookie);
				add(glow_server + "ping_glow");
			}
		};
		if (runningOnLinux) {
			try {
				String response = "";
				Process process = new ProcessBuilder(command).start();
				BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(process.getInputStream()));
				String line = null;
				while ((line = bufferedReader.readLine()) != null) {
					response += line;
				}
				process.waitFor();
				return response.contains("Ping saved");
			} catch (IOException | InterruptedException e) {
				e.printStackTrace();
			}
			return false;
		} else {
			Shell ssh = establishSsh();
			try {
				String response = new Shell.Plain(ssh).exec(convertCommandToString(command));
				return response.contains("Ping saved");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return false;
		}
	}

	protected boolean renewCookieExpiration() {
		/*
		 * If the ping doesn't work, tries to renew cookie
		 */
		if (!pingGlow()) {
			requestCookie();
		} else {
			// If the ping worked, just returns true
			return true;
		}
		// If the cookie was attempted to be renewed, tries to ping again and
		// returns the result for the caller.
		return pingGlow();
	}

	/*
	 * FUTURE private void upsertExperiment() { // TO DO }
	 */

	public void logout() {
		removeCookie();
		username = "";
		password = "";
	}

	private void removeCookie() {
		if (runningOnLinux) {
			File c = new File(this.cookie);
			try {
				Files.delete(c.toPath());
			} catch (Exception e) {
				System.out.println("Error: Unable to delete program cookie");
			}
		} else {
			Shell ssh = establishSsh();
			try {
				try {
					new Shell.Plain(ssh).exec("rm " + cookie);
				} catch (NullPointerException e) {
					// Meh no cookie getting removed because they logged out or
					// something. Cookie should already be removed if they
					// logged out too, so this would only trigger on shutdown.
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	@Override
	protected Object call() throws Exception {
		boolean serverUp = pingGlow();
		if (serverUp) {
			connected = true;
		} else {
			System.out.println("Error contacting the server.");
		}
		return null;
	}
}
