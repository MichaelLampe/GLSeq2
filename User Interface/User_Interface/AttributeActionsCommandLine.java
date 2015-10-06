package application;

import java.util.ArrayList;

public final class AttributeActionsCommandLine extends AttributeActions {

	public void setAttributes(String[] inputArgs) {
		ArrayList<Attribute> attributeArray = new ArrayList<Attribute>(
				Attributes.getInstance().attributesCollection.values());
		for (String param : inputArgs) {
			String[] parts = param.split("=");
			// Complexity is fine when it is cheap.
			for (Attribute att : attributeArray) {
				if (att.getUiName().equals(parts[0])) {
					Attribute current = Attributes.getInstance().attributesCollection.get(att.getName());
					try {
						current.setValue(parts[1]);
					} catch (Exception e) {
						// Improper field was entered
					}
				}
			}
		}
	}

	public void returnJson() {
		System.out.println("{\"attributes\":[");
		ArrayList<Attribute> attributeArray = new ArrayList<Attribute>(
				Attributes.getInstance().attributesCollection.values());
		for (int i = attributeArray.size() - 1; i >= 0; i--) {
			if (attributeArray.get(i).getCategory().equals(Category.SCRIPT.name)) {
				attributeArray.remove(i);
			}
		}
		for (int i = 0; i < attributeArray.size(); i++) {
			System.out.print(formatJson(attributeArray.get(i).getUiName()));
			if (i < attributeArray.size() - 1) {
				System.out.print(",");
			}
			System.out.println("");
		}
		System.out.println("]}");
	}

	/*
	 * JSON formatted String to return to GLOW
	 */
	private String formatJson(String field_name) {
		String jsonKey = "{\"key\":";
		AttributesJSON current = Enum.valueOf(AttributesJSON.class, field_name);
		jsonKey += "\"" + current.name() + "\"";
		if (current.category != null) {
			jsonKey += ", \"category\":\"" + current.category + "\"";
		}
		if (current.options != null) {
			jsonKey += ", \"options\":[" + current.options + "]";
		}
		if (current.defaultVal != null) {
			jsonKey += ", \"default\":\"" + defaultNullIsEmpty(current.defaultVal) + "\"";
		}
		if (current.description != null) {
			jsonKey += ", \"description\":\"" + current.description + "\"";
		}
		jsonKey += "}";
		return jsonKey;
	}

	private String defaultNullIsEmpty(String defaultVal) {
		if (defaultVal.equals("NULL")) {
			return "";
		}
		return defaultVal;
	}
}
