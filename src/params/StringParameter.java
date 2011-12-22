package params;

public class StringParameter extends Parameter {
	String value = null;
	public StringParameter(String key, String name, String description) {
		super(key, name, description);
	}

	public StringParameter defaultValue(String defaultValue)
	{
		this.value = defaultValue;
		super.setOptional();
		return this;
	}
	
	@Override
	public String parse(String value) {
		this.value = value.trim();
		return null;
	}

	@Override
	public String getValueAsString() {
		return (value == null ? "null" : value);
	}

}
