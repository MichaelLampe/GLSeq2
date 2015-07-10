package Tests;

import org.glbrc.glseq2.Project_Files.AttributesJSON;
import org.glbrc.glseq2.Project_Files.Attributes;

public class JsonTest {
  {

  }

  static boolean runTest() {
    Attributes att = new Attributes();
    // The -1 is here because we don't present the RunID to the user, even
    // though we do want to embed that into the Attribute file if possible.
    if (AttributesJSON.size == att.getFieldNumber() - 1) {
      System.out.println("JsonTests passed");
      return true;
    } else {
      System.out.println("JsonTests failed");
      System.out.println("Attributes JSON Size:" + AttributesJSON.size);
      System.out.println("Attributes Class Field Number:" + att.getFieldNumber());
      return false;
    }
  }
}
