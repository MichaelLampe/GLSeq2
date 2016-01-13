package application;

import java.util.ArrayList;

public final class AttributeActionsCommandLine extends AttributeActions {

	public void setAttributes(String[] inputArgs) {
		ArrayList<AttributeFactorySingleton.Attribute> attributeArray = AttributeFactorySingleton
				.getInstance().getAllAttributes();
		for (String param : inputArgs) {
			String[] parts = param.split("=");
			// Complexity is fine when it is cheap.
			for (AttributeFactorySingleton.Attribute att : attributeArray) {
				if (att.getUiName().equals(parts[0])) {
					try {
						AttributeFactorySingleton.getInstance()
								.setAttributeValue(att.getName(), parts[1]);
					} catch (NoSuchKeyInAttributeFactoryException e1) {
						/*
						 * Did not find attribute key.
						 */
					}
				}
			}
		}
	}

	public void returnJson() {
		System.out.println("{\"attributes\":[");
		ArrayList<AttributeFactorySingleton.Attribute> attributeArray = AttributeFactorySingleton
				.getInstance().getAllAttributes();
		for (int i = attributeArray.size() - 1; i >= 0; i--) {
			/*
			 * We want to give them all the categories that are not in the
			 * script.
			 */
			if (attributeArray.get(i).getCategory().equals("script")) {
				attributeArray.remove(i);
			}
		}
		for (int i = 0; i < attributeArray.size(); i++) {
			/*
			 * They have SSH'd into this file, so they can read the system out
			 * stream.
			 */
			try {
				System.out.print(formatJson(attributeArray.get(i).getUiName()));
			} catch (NoSuchKeyInAttributeFactoryException e) {
				/*
				 * Key wasn't there so don't print it.
				 */
			}
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
	private String formatJson(String uiName)
			throws NoSuchKeyInAttributeFactoryException {
		String jsonKey = "{\"key\":";
		AttributeFactorySingleton.Attribute current;
		current = AttributeFactorySingleton.getInstance().getAttribute(uiName);
		jsonKey += "\"" + current.getUiName() + "\"";
		if (current.getCategory() != null) {
			jsonKey += ", \"category\":\"" + current.getCategory() + "\"";
		}
		if (current.getOptions() != null) {
			jsonKey += ", \"options\":[" + current.getOptionsPresplit() + "]";
		}
		if (current.getDefault() != null) {
			jsonKey += ", \"default\":\""
					+ defaultNullIsEmpty(current.getDefault()) + "\"";
		}
		if (current.getToolTip() != null) {
			jsonKey += ", \"description\":\"" + current.getToolTip() + "\"";
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
