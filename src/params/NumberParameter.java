package params;

public abstract class NumberParameter<T extends Number> extends Parameter {
	protected T minValue; // inclusive
	protected T maxValue;	  // exclusive
	protected T value;
	
	public NumberParameter(String key, String name, String description) {
		super(key, name, description);
	}

	public NumberParameter<T> defaultValue(T defaultValue)
	{
		value = defaultValue;
		super.setOptional();
		return this;
	}
	
	public NumberParameter<T> minValue(T minValue) 
	{
		this.minValue = minValue;
		return this;
	}

	public NumberParameter<T> maxValue(T maxValue) 
	{
		this.maxValue = maxValue;
		return this;
	}
	
	@Override
	public abstract String parse(String value);

	@Override
	public String getValueAsString() {
		return String.valueOf(value);
	}

	public T getValue()
	{
		return value;
	}
}
