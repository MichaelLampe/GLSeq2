package application;

class Attribute implements Attributable{

  private final String name;
  private final String uiname;
  private final String toolTip;
  private String value; 
  
  Attribute(String name,String uiname, String toolTip){
    this.uiname = uiname;
    this.name = name;
    this.toolTip = toolTip;
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
