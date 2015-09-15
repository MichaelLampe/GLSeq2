package application;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.util.ArrayList;

import javafx.concurrent.Task;

public class GlowRequest extends Task<Object> {
  protected final String glow_server = "https://glow-trunk.glbrc.org/";

  protected String username;
  // While a string is not the best way to store this, the password
  // Will be logged in any string we make that calls a curl command anyway.
  protected String password = "";
  protected String cookie;
  // If connected to GLOW by cookie and ping was true
  protected boolean connected = false;

  public GlowRequest(String username, String password) {
    this.username = username;
    this.password = password;
    this.cookie = "cjar_" + username;
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

  boolean requestCookie() {

    @SuppressWarnings("serial")
    ArrayList<String> command = new ArrayList<String>() {
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
      Process process = new ProcessBuilder(command).directory(
          new File(System.getProperty("user.dir"))).start();
      BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(
          process.getInputStream()));
      String line = null;
      while ((line = bufferedReader.readLine()) != null) {
        response += line;
      }
      process.waitFor();
      // The actual response if correct should be about 100, but this future
      // proofs it a bit.
      // If the password is wrong it returns the whole html page so even 200
      // should be pretty safe.
      // The greater than protects us against incorrectly directed requests that
      // don't do anything
      // Nor get anything back making the program think they are correct
      return (response.length() < 200 && response.length() > 2);
    } catch (IOException | InterruptedException e) {
      // Don't print stack trace here as it would reveal the user's password in
      // that trace
    }
    return false;
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
    try {
      String response = "";
      Process process = new ProcessBuilder(command).start();
      BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(
          process.getInputStream()));
      String line = null;
      while ((line = bufferedReader.readLine()) != null) {
        response += line;
        System.out.println(line);
      }
      process.waitFor();
      return response.contains("Ping saved");
    } catch (IOException | InterruptedException e) {
      e.printStackTrace();
    }
    return false;
  }

  protected boolean renewCookieExpiration() {
    // If the ping doesn't work, tries to renew cookie
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

  private void upsertExperiment() {
    // TO DO
  }

  private void runCommand(String command) {
    // Process builder
  }

  public void logout() {
    removeCookie();
    username = "";
    password = "";
  }

  private void removeCookie() {
    File c = new File(this.cookie);
    try {
      Files.delete(c.toPath());
    } catch (Exception e) {
      System.out.println("Error: Unable to delete program cookie");
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
