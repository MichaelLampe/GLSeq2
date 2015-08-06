package application;

import java.util.LinkedHashMap;
import java.util.Map;
public class Attributes {
  
  protected Map<String,Attribute> attributesCollection;
  
  private static Attributes instance = null;
  
  protected Attributes(){
    attributesCollection = new LinkedHashMap<String,Attribute>();
    setupFields();
  }  
  public static Attributes getInstance(){
    if(instance == null){
      instance = new Attributes();
    }
    return instance;
  }
  /*
   * Initializes all the fields with the name and default val
   * that they have in the enum.
   */
  private void setupFields(){
    for (AttributesJSON field : AttributesJSON.values()){
      Attribute currentAttribute = new Attribute(field.field_name,field.name(),field.description,field.category);
      currentAttribute.setValue(field.defaultVal);
      attributesCollection.put(currentAttribute.getName(), currentAttribute);
    }
  }
}
