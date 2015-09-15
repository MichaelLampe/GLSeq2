package application;

class Attribute implements Attributable{

  private final String name;
  private final String uiname;
  private final String toolTip;
  private final String category;
  private final String defaultValue;
  private String value; 
  
  Attribute(String name,String uiname, String toolTip,String category,String defaultValue){
    this.name = name;
    this.uiname = uiname;
    this.toolTip = toolTip;
    this.category = category;
    this.defaultValue = defaultValue;
  }
  public String getCategory(){
    return category;
  }
  @Override
  public String getValue() {
    return value;
  }
  @Override
  public String getToolTip(){
    return toolTip;
  }
  @Override
  public String getName(){
    return name;
  }
  @Override
  public void setValue(String propertyUpdate) {
    this.value = propertyUpdate;
  }
  public String getDefault(){
    return defaultValue;
  }
  public String getUiName(){
    return uiname;
  }
  
  @Override
  public String toString(){
    String message = "";
    message += "Name:" + name + " with a value of " + value;
    return message;  
  }
}
