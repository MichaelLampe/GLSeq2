package application;

import java.io.IOException;
import javafx.application.Application;
import javafx.fxml.FXMLLoader;
import javafx.scene.Scene;
import javafx.scene.layout.VBox;
import javafx.stage.Stage;

public class Main extends Application {

  public static Scene scene;
  private static final AttributeActions att = new AttributeActions();

  @Override
  public void start(Stage primaryStage) {
    try {
      VBox root = (VBox) FXMLLoader.load(Main.class.getResource("MainPage.fxml"));
      String cssFile = this.getClass().getResource("MainPage.css").toExternalForm();
      scene = new Scene(root, 873, 770);
      scene.getStylesheets().add(cssFile);
      primaryStage.setResizable(false);
      primaryStage.setScene(scene);
      primaryStage.setTitle("GLSeq2 User Interface");
      primaryStage.show();
      // Update attributes from previous save
      att.setAttributes();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {
    // Command line args were passed
    if (args.length > 0) {
      // Returns a JSON format of all the attribute file values
      if (args[0].equals("JSON")) {
        returnJson();
        // Just makes an attribute file.
      } else {
        generateAttributeFile(args);
      }
      // Launch the UI as normal
    } else {
      launch(args);
    }
  }

  /*
   * Prints out a JSON formatted version of all the attribute file fields to
   * STDIO (Which can be captured by the caller).
   */
  private static void returnJson() {
    att.returnJson();
    System.exit(0);
  }

  /*
   * Create a new attribute file based on the args. FILE_NAME= changes the
   * attribute file name FILE_LOCATION= changes where the file is saved All
   * other args are field names, = , and then the value that should be placed in
   * that field. Note: There is currently very little validation that goes on
   * here, as all the params are strings. May look into this in the future.
   */
  private static void generateAttributeFile(String[] args) {
    try {
      String file_name = null;
      String file_location = null;
      att.setAttributes(args);
      try {
        if (args[0].contains("FILE_NAME")) {
          file_name = args[0].split("=")[1];
          if (args[1].contains("FILE_LOCATION")) {
            file_location = args[1].split("=")[1];
          }
        } else {
          if (args[0].contains("FILE_LOCATION")) {
            file_location = args[0].split("=")[1];
          }
        }
      } catch (IndexOutOfBoundsException e) {
        // Do nothing here.
      }
      System.out
          .println("ATTIRUBTE_FILE_PATH=" + att.writeAttributesFile(file_name, file_location));
    } catch (IOException e) {
      System.out
          .println("Error constructing attribute file. Likely unable to create attribute file where the JAR file exists.");
    }
    // Exit program
    System.exit(0);
  }
}