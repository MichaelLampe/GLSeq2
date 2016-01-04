package application;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import javafx.application.Application;
import javafx.fxml.FXMLLoader;
import javafx.scene.Scene;
import javafx.scene.layout.VBox;
import javafx.stage.Stage;

public class Main extends Application {

	public static Scene scene;
	public static Stage stage;
	private static final AttributeActionsUI att = new AttributeActionsUI();

	@Override
	public void start(Stage primaryStage) throws IOException {
		VBox root = (VBox) FXMLLoader.load(Main.class.getResource("/application/MainPage.fxml"));
		String cssFile = this.getClass().getResource("/application/MainPage.css").toExternalForm();
		scene = new Scene(root, 873, 770);
		scene.getStylesheets().add(cssFile);
		primaryStage.setResizable(false);
		primaryStage.setScene(scene);
		primaryStage.setTitle("|   GLSeq2 User Interface   -   Not logged into GLOW   |");
		stage = primaryStage;
		primaryStage.show();
		// Update attributes from previous save
		String defaultAttributeFile = "AttributesConfig.txt";
		try {
			att.setAttributesTextConfig(new File(defaultAttributeFile));
		} catch (FileNotFoundException e) {
			System.out.println("Presaved attributes config file not found (" + defaultAttributeFile + ").");
		}
	}

	public static void main(String[] args) {
		launch(args);
	}
}