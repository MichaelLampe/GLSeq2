package application;

public enum Category {
	DATA("Data and Library"), REFERENCE("Reference Attributes"), RUN("Run Attributes"), PREPROCESS(
			"Pre-Processing Attributes"), COMMON_PROCESS("Common Processing Options"), RSEM(
					"RSEM Attributes"), ENVIRONMENT("Environment Attributes"), SCRIPT("Script Args");
	public final String name;

	Category(String name) {
		this.name = name;
	}
}