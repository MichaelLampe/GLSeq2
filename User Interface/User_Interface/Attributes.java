package application;

import java.util.LinkedHashMap;
import java.util.Map;

public class Attributes {

	protected Map<String, Attribute> attributesCollection;

	private static Attributes instance = null;

	protected Attributes() {
		attributesCollection = new LinkedHashMap<String, Attribute>();
		setupFields();
	}

	public static Attributes getInstance() {
		if (instance == null)
			instance = new Attributes();
		return instance;
	}

	/*
	 * Initializes all the fields with the name and default val that they have
	 * in the enum.
	 */
	private void setupFields() {
		// This is a special field we will only allow in the UI
		Attribute alignmentSpecialOptions = new Attribute(
				"alignmentSpecialOptions", "alignmentSpecialOptions", "", "",
				"", "");
		attributesCollection.put(alignmentSpecialOptions.getName(), alignmentSpecialOptions);
		
		Attribute htSeqSpecialOptions = new Attribute("HtSeqSpecialOptions",
				"HtSeqSpecialOptions", "", "", "", "");
		attributesCollection.put(htSeqSpecialOptions.getName(), htSeqSpecialOptions);
		
		Attribute rsemSpecialOptions = new Attribute("RsemSpecialOptions",
				"RsemSpecialOptions", "", "", "", "");
		attributesCollection.put(rsemSpecialOptions.getName(), rsemSpecialOptions);
		
		Attribute featureCountsSpecialOptions = new Attribute(
				"FeatureCountsSpecialOptions", "FeatureCountsSpecialOptions",
				"", "", "", "");
		attributesCollection.put(featureCountsSpecialOptions.getName(), featureCountsSpecialOptions);
		
		Attribute cufflinksSpecialOptions = new Attribute(
				"CufflinksSpecialOptions", "CufflinksSpecialOptions", "", "",
				"", "");
		attributesCollection.put(cufflinksSpecialOptions.getName(), cufflinksSpecialOptions);

		
		for (AttributesJSON field : AttributesJSON.values()) {
			Attribute currentAttribute = new Attribute(field.field_name,
					field.name(), field.description, field.category,
					field.defaultVal, field.options);
			currentAttribute.setValue(field.defaultVal);
			attributesCollection.put(currentAttribute.getName(),
					currentAttribute);
		}
	}
}
