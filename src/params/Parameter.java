package params;

public abstract class Parameter {
	private String key;
	private String name;
	private String description;
	private boolean isOptional;
	
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
}
