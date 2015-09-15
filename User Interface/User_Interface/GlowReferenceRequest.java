package application;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import javafx.application.Platform;

public class GlowReferenceRequest extends GlowRequest {
  private String reference_id;

  public GlowReferenceRequest(GlowRequest connection, String reference_id) {
    super(connection.getUsername(), connection.getPassword());
    this.reference_id = reference_id;
    this.connected = connection.getConnected();
  }

  private List<String> requestReferenceFiles(String reference_id) {
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
    if (renewCookieExpiration()) {
      Process processFiles = new ProcessBuilder(requestReferenceFiles(reference_id)).start();
      Process processName = new ProcessBuilder(requestReferenceName(reference_id)).start();

      BufferedReader fileR = new BufferedReader(
          new InputStreamReader(processFiles.getInputStream()));
      BufferedReader nameR = new BufferedReader(new InputStreamReader(processName.getInputStream()));

      // Parse for GFF and Fasta file lines
      List<String> commands = new ArrayList<String>();
      commands
          .addAll(fileR
              .lines()
              .filter(
                  line -> (line.contains("<td>GFF</td>"))
                      || (line.contains("<td>Sequence:FastA</td>"))).collect(Collectors.toList()));

      commands.addAll(nameR.lines().map(line -> line.split(","))
          // Flat map makes it happy that we split it again.
          .flatMap(Arrays::stream).filter(line -> (line.contains("name")))
          .collect(Collectors.toList()));
      // Parse for reference name
      processFiles.waitFor();
      processName.waitFor();
      for (String line : commands) {
        if (line.contains("name")) {
          String[] temp = line.split("\"");
          String refname = temp[3];
          Attributes.getInstance().attributesCollection.get(
              AttributesJSON.referenceGenome.field_name).setValue(refname);
        } else if (line.contains("<td>GFF</td>")) {
          String[] temp = line.split("class=\"detail_link\">");
          String gff = temp[1].split("</a>")[0];
          Attributes.getInstance().attributesCollection.get(AttributesJSON.referenceGff.field_name)
              .setValue(gff);
        } else {
          String[] temp = line.split("class=\"detail_link\">");
          String fasta = temp[1].split("</a>")[0];
          Attributes.getInstance().attributesCollection.get(
              AttributesJSON.referenceFasta.field_name).setValue(fasta);
        }
      }
      // Only javafx thread can update javafx
      Platform.runLater(() -> { 
        UpdateUI.getInstance().updateDefaults();
      });
    } else {
      System.out.println("Did not run because no connection to glow was established.");
    }
    return null;
  }
}
