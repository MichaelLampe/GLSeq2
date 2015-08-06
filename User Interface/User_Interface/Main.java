package application;

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
      // ok
    }
  }

  public static void main(String[] args) {
    launch(args);
  }
}