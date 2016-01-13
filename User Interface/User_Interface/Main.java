package application;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import javafx.application.Application;
import javafx.fxml.FXMLLoader;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.layout.VBox;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.stage.Stage;

public class Main extends Application {

	private static final String INITIAL_TITLE = "|   GLSeq2 User Interface   -   Not logged into GLOW   |";
	private static final String ATTRIBUTES_CONFIG_FILE = "AttributesConfig.txt";
	private static Scene scene;
	private static Stage stage;
	private static final AttributeConfigurationFile att = new AttributeConfigurationFile();

	/*
	 * Sets the stage title
	 */
	public static void setStageTitle(String title) {
		stage.setTitle(title);
	}

	/*
	 * Returns stage title back to INITIAL_TITLE
	 */
	public static void setStageTitleDefault() {
		stage.setTitle(INITIAL_TITLE);
	}

	/*
	 * Opens a file chooser dialog and returns the chosen file.
	 */
	public static File openFileChooser(FileChooser fileChooser) {
		return fileChooser.showOpenDialog(stage);
	}

	/*
	 * Opens a directory chooser dialog and returns the chosen directory.
	 */
	public static File openDirectoryChooser(DirectoryChooser directoryChooser) {
		return directoryChooser.showDialog(stage);
	}

	public static Node lookupObjectById(String id) throws SceneNotYetReadyException {
		if (scene == null) {
			throw new SceneNotYetReadyException();
		}
		return scene.lookup(id);
	}

	@Override
	public void start(Stage primaryStage) throws IOException {
		/*
		 * Basic class loading stuff and setting application size.
		 */
		VBox root = (VBox) FXMLLoader.load(Main.class
				.getResource("/application/MainPage.fxml"));
		String cssFile = this.getClass()
				.getResource("/application/MainPage.css").toExternalForm();
		scene = new Scene(root, 873, 770);

		/*
		 * Scene is assigned, update defaults.
		 */
		UpdateUserInterfaceSingleton.getInstance().updateDefaults();
		scene.getStylesheets().add(cssFile);
		/*
		 * It doesn't scale relative to the resize, so it doesn't look good when
		 * it is able to be resized.
		 */
		primaryStage.setResizable(false);
		primaryStage.setScene(scene);
		primaryStage.setTitle(INITIAL_TITLE);
		stage = primaryStage;
		primaryStage.show();

		try {
			att.setAttributesTextConfig(new File(ATTRIBUTES_CONFIG_FILE));
		} catch (FileNotFoundException e) {
			// Just move on.
		}
	}

	public static void main(String[] args) {
		// Launch UI
		launch(args);
	}
}