package params;

public class StringParameter extends Parameter {
	String value;
	protected StringParameter(String key, String name, String description) {
		super(key, name, description);
	}

	public StringParameter defaultValue(String defaultValue)
	{
		this.value = defaultValue;
		super.setOptional();
		return this;
	}
	
	@Override
	public boolean parse(String value) {
		return false;
	}

	@Override
	public String getValueAsString() {
		// TODO Auto-generated method stub
		return null;
	}

}
