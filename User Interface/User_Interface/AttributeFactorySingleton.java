package application;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.xml.sax.SAXException;

public class AttributeFactorySingleton {

	/*
	 * Holds all the attributes
	 */
	private Map<String, Attribute> attributesCollection;
	/*
	 * Contains the attributes
	 */
	
	private static final String ATTRIBUTES_XML_ROOT_NAME = "Attributes";
	private static final String ATTRIBUTES_XML = "xml_data/Attributes.xml";
	private static final String ATTRIBUTES_XML_FIELD_NAME_ATTRIBUTE = "name";
	private static final String ATTRIBUTES_XML_CATEGORY_ATTRIBUTE = "category";
	private static final String ATTRIBUTES_XML_OPTIONS_ATTRIBUTE = "options";
	private static final String ATTRIBUTES_XML_DEFAULT_VALUE_ATTRIBUTE = "default_value";
	private static final String ATTRIBUTES_XML_DESCRIPTION_ATTRIBUTE = "description";
	
	public static final String ALIGNMENT_SPECIAL_OPTIONS_NAME = "alignmentSpecialOptions";
	public static final String HTSEQ_SPECIAL_OPTIONS_NAME ="HtSeqSpecialOptions";
	public static final String RSEM_SPECIAL_OPTIONS_NAME ="RsemSpecialOptions";
	public static final String FEATURE_COUNTS_SPECIAL_OPTIONS_NAME ="FeatureCountsSpecialOptions";
	public static final String CUFFLINKS_SPECIAL_OPTIONS_NAME ="CufflinksSpecialOptions";
	
	private static AttributeFactorySingleton instance = null;

	protected AttributeFactorySingleton() {
		attributesCollection = new LinkedHashMap<String, Attribute>();

		/*
		 * This is used to setup the special option fields.
		 */
		setupSpecialOptionFields();

		/*
		 * Setup the attributes based on the Attributes.xml file.
		 */
		try {
			createAttributes(ATTRIBUTES_XML);
		} catch (ParserConfigurationException | SAXException | IOException e) {
			System.out.println("Unable to create attributes from XML file.");
		}
	}

	/*
	 * Holds all the attributes, so it is a singleton providing one access
	 * point.
	 */
	public static AttributeFactorySingleton getInstance() {
		if (instance == null)
			instance = new AttributeFactorySingleton();
		return instance;
	}

	public Attribute getAttribute(String attributeName) throws NoSuchKeyInAttributeFactoryException {
		if (containsAttributeKey(attributeName)) {
			return attributesCollection.get(attributeName);
		} else{
			throw new NoSuchKeyInAttributeFactoryException();
		}
	}

	public ArrayList<Attribute> getAllAttributes() {
		return new ArrayList<Attribute>(attributesCollection.values());
	}

	public void setAttributeValue(String attributeName, String value) throws NoSuchKeyInAttributeFactoryException {
		Attribute attributeToChange = getAttribute(attributeName);
		attributeToChange.setValue(value);
	}

	public boolean containsAttributeKey(String keyName) {
		return attributesCollection.containsKey(keyName);
	}

	/*
	 * Initializes all the fields with the name and default val that they have
	 * in the enum.
	 */
	private void setupSpecialOptionFields() {
		// This is a special field we will only allow in the UI
		Attribute alignmentSpecialOptions = new Attribute(
				ALIGNMENT_SPECIAL_OPTIONS_NAME, ALIGNMENT_SPECIAL_OPTIONS_NAME, "", "",
				"", "");
		attributesCollection.put(alignmentSpecialOptions.getName(),
				alignmentSpecialOptions);

		Attribute htSeqSpecialOptions = new Attribute(HTSEQ_SPECIAL_OPTIONS_NAME,
				HTSEQ_SPECIAL_OPTIONS_NAME, "", "", "", "");
		attributesCollection.put(htSeqSpecialOptions.getName(),
				htSeqSpecialOptions);

		Attribute rsemSpecialOptions = new Attribute(RSEM_SPECIAL_OPTIONS_NAME,
				RSEM_SPECIAL_OPTIONS_NAME, "", "", "", "");
		attributesCollection.put(rsemSpecialOptions.getName(),
				rsemSpecialOptions);

		Attribute featureCountsSpecialOptions = new Attribute(
				FEATURE_COUNTS_SPECIAL_OPTIONS_NAME, FEATURE_COUNTS_SPECIAL_OPTIONS_NAME,
				"", "", "", "");
		attributesCollection.put(featureCountsSpecialOptions.getName(),
				featureCountsSpecialOptions);

		Attribute cufflinksSpecialOptions = new Attribute(
				CUFFLINKS_SPECIAL_OPTIONS_NAME, CUFFLINKS_SPECIAL_OPTIONS_NAME, "", "",
				"", "");
		attributesCollection.put(cufflinksSpecialOptions.getName(),
				cufflinksSpecialOptions);
	}

	private void createAttributes(String xmlPath)
			throws ParserConfigurationException, SAXException, IOException {
		/*
		 * Load XML file in via xmlPath
		 */
		DocumentBuilder xml = DocumentBuilderFactory.newInstance()
				.newDocumentBuilder();
		Document attributes;
		try {
			attributes = xml.parse(new File(xmlPath));
		} catch (FileNotFoundException e) {
			attributes = xml.parse(TreeViewOptionsLoader.class
					.getResourceAsStream(xmlPath));
		}

		/*
		 * Go through the XML document loaded in via xmlPath
		 */
		Node attributeNode = attributes.getElementsByTagName(ATTRIBUTES_XML_ROOT_NAME)
				.item(0);
		for (int i = 0; i < attributeNode.getChildNodes().getLength(); i++) {

			// Grab the current xml element
			Node xmlAttribute = attributeNode.getChildNodes().item(i);

			/*
			 * #text nodes are not useful to us, so we exclude them and just
			 * process non #text named nodes.
			 */
			if (!xmlAttribute.getNodeName().equals("#text")) {
				String name = xmlAttribute.getNodeName();
				String field_name = "";
				String category = "";
				String options = "";
				String default_value = "";
				String description = "";

				try {
					field_name = xmlAttribute.getAttributes()
							.getNamedItem(ATTRIBUTES_XML_FIELD_NAME_ATTRIBUTE).getNodeValue();
				} catch (NullPointerException e) {
					// No field name value in XML node
				}

				try {
					category = xmlAttribute.getAttributes()
							.getNamedItem(ATTRIBUTES_XML_CATEGORY_ATTRIBUTE).getNodeValue();
				} catch (NullPointerException e) {
					// No category value in XML node
				}

				try {
					options = xmlAttribute.getAttributes()
							.getNamedItem(ATTRIBUTES_XML_OPTIONS_ATTRIBUTE).getNodeValue();
				} catch (NullPointerException e) {
					// No options value in XML node
				}

				try {
					default_value = xmlAttribute.getAttributes()
							.getNamedItem(ATTRIBUTES_XML_DEFAULT_VALUE_ATTRIBUTE).getNodeValue();
				} catch (NullPointerException e) {
					// No default value in XML node
				}

				try {
					description = xmlAttribute.getAttributes()
							.getNamedItem(ATTRIBUTES_XML_DESCRIPTION_ATTRIBUTE).getNodeValue();
				} catch (NullPointerException e) {
					// No description value in XML node
				}
				Attribute currentAttribute = new Attribute(field_name, name,
						description, category, default_value, options);
				currentAttribute.setValue(default_value);

				/*
				 * Hash the attribute to the attribute's UI name for later
				 * lookup. This is also the name of the node from the XML
				 * document. There are a few hard-coded UI names (By switching
				 * from enum -> XML they had to be hard-coded), so UI names
				 * should not be changed unless you want to try to find those.
				 */
				attributesCollection.put(name, currentAttribute);
			}
		}
	}

	/*
	 * AttributeFactory defines what an attribute is
	 */
	class Attribute {
		private final String name;
		private final String uiname;
		private final String toolTip;
		private final String category;
		private final String defaultValue;
		private String value;
		private final String optionsPresplit;
		private final String[] options;

		Attribute(String name, String uiname, String toolTip, String category,
				String defaultValue, String options) {
			this.name = name;
			this.uiname = uiname;
			this.toolTip = toolTip;
			this.category = category;
			this.defaultValue = defaultValue;

			// Parse options on being added
			if (options != null) {
				this.options = options.split(",");
			} else {
				this.options = null;
			}

			/*
			 * Used by GLOW to get all options. They do the parsing on their
			 * end.
			 */
			this.optionsPresplit = options;
		}

		public String[] getOptions() {
			return this.options;
		}

		public String getCategory() {
			return this.category;
		}

		public String getValue() {
			return this.value;
		}

		public String getToolTip() {
			return this.toolTip;
		}

		public String getName() {
			return this.name;
		}

		public void setValue(String propertyUpdate) {
			this.value = propertyUpdate;
		}

		public String getDefault() {
			return this.defaultValue;
		}

		public String getUiName() {
			return this.uiname;
		}

		public String getOptionsPresplit() {
			return this.optionsPresplit;
		}
	}
}
