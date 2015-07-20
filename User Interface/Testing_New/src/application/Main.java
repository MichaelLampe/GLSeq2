package application;
	
import javafx.application.Application;
import javafx.fxml.FXMLLoader;
import javafx.scene.Scene;
import javafx.scene.layout.VBox;
import javafx.stage.Stage;

public class Main extends Application {
	@Override
	public void start(Stage primaryStage) {
		try {
		    VBox root = (VBox) FXMLLoader.load(Main.class.getResource("MainPage.fxml"));
			Scene scene = new Scene(root,873,770);
			primaryStage.setResizable(false);
			primaryStage.setScene(scene);
			primaryStage.setTitle("GLSeq2 User Interface");
			primaryStage.show();
		} catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		launch(args);
	}
}
