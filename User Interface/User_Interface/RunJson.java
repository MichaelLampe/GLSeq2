package application;

public class RunJson {
  private static final AttributeActions att = new AttributeActions();

  // Can be launched without the other parts of the program by using
  // -cp as the java option and referring to this as application/RunJson
  public static void main(String[] args) {
    returnJson();
    System.exit(0);
  }

  /*
   * Prints out a JSON formatted version of all the attribute file fields to
   * STDIO (Which can be captured by the caller).
   */
  private static void returnJson() {
    att.returnJson();
  }
}
