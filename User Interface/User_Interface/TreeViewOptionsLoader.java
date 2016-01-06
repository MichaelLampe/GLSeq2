package application;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;

import javafx.scene.control.CheckBoxTreeItem;
import javafx.scene.control.TreeItem;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.xml.sax.SAXException;

public class TreeViewOptionsLoader {

	private Document myOptions;
	private static HashMap<String,String> commandMap = new HashMap<String, String>();
	
	public TreeViewOptionsLoader() {
		DocumentBuilder xml = null;
		try {
			xml = DocumentBuilderFactory.newInstance().newDocumentBuilder();
		} catch (ParserConfigurationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		try {
			String xml_path = "xml_data/treeOptions.xml";
			// This should fix running in eclipse vs JAR file loading of the XML file
			try {
				myOptions = xml.parse(new File(xml_path));
			} catch (FileNotFoundException e) {
				myOptions = xml.parse(TreeViewOptionsLoader.class
						.getResourceAsStream(xml_path));
			}
		} catch (SAXException | IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public TreeItem<String> createTreeElementfromXml(String itemName)
			throws DuplicateElementsInXmlError {
		TreeItem<String> optionsHead = new TreeItem<String>("Options");

		// This is here because we could get some hidden bugs if we don't notice
		// that we are grabbing multiple XML elements, invalidating a key
		// assumption that this code has.
		if (myOptions.getElementsByTagName(itemName).getLength() > 1) {
			throw new DuplicateElementsInXmlError();
		}
		Node itemNode = myOptions.getElementsByTagName(itemName).item(0);
		Node optionsNode = itemNode.getChildNodes().item(1);

		// First and last nodes are just text, so we ignore those two.
		for (int i = 1; i < optionsNode.getChildNodes().getLength() - 1; i++) {
			if (!optionsNode.getChildNodes().item(i).getNodeName()
					.equals("#text")) {
				optionsHead.getChildren().add(
						parseSingleOption(itemName, optionsNode.getChildNodes().item(i)));
			}
		}
		return optionsHead;
	}

	private TreeItem<String> parseSingleOption(String itemName, Node optionXml) {
		// Every option has three properties:

		// Command
		// Input Type
		// Default Value
		String name = itemName + "_" + optionXml.getNodeName();
		
		CheckBoxTreeItem<String> option = new CheckBoxTreeItem<String>(
				name);

		String input_type = optionXml.getAttributes()
				.getNamedItem("input_type").getNodeValue();
		String default_value = optionXml.getAttributes()
				.getNamedItem("default_value").getNodeValue();
		

		
		// We'll use this later as we don't necessarily want to display this.
		getCommandMap().put(name, optionXml.getAttributes().getNamedItem("command").getNodeValue());
		if (!input_type.equals("add")) {
			option.getChildren()
					.add(new SpecialInputTreeItem<String>(input_type,
							default_value));
		}
		return option;
	}

	public static HashMap<String,String> getCommandMap() {
		return commandMap;
	}

	public static void setCommandMap(HashMap<String,String> commandMap) {
		TreeViewOptionsLoader.commandMap = commandMap;
	}

}
