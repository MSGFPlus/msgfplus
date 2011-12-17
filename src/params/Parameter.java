package params;

public abstract class Parameter {
	private String key;
	private String name;
	private String description;
	private boolean isOptional;
	private boolean isValueAssigned = false;
	
	protected Parameter(String key, String name, String description, boolean isOptional)
	{
		this.key = key;
		this.name = name;
		this.description = description;
		this.isOptional = isOptional;
	}
	
	public String getKey()	{ return key; }
	public String getName()			{ return name; }
	public String getDescription()	{ return description; }
	public boolean isOptional()		{ return isOptional; }
	
	public String toString()
	{
		String usage = "-" + key + " " + name;
		if(isOptional)
			usage = "[" + usage + "]";
		usage = usage + " " + "(" + description + ")";
		return usage;
	}
	
	public void setValueAssigned() { this.isValueAssigned = true; }
	public boolean isValueAssigned() { return isValueAssigned; }
	public boolean isValid() 
	{ 
		if(!isOptional && !isValueAssigned())
			return false;
		return true; 
	}
	
	public abstract boolean parse(String value);
	public abstract String getValueAsString();
}
