package edu.ucsd.msjava.params;

public abstract class RangeParameter<T extends Comparable<T>> extends Parameter {
	protected T min = null;	
	protected T max = null;
	protected T minValue; 		// default: inclusive
	protected T maxValue;		// default: exclusive
	protected boolean isMinInclusive = true;
	protected boolean isMaxInclusive = false;
	
	public RangeParameter(String key, String name, String description) {
		super(key, name, description);
	}

	public RangeParameter<T> minValue(T minValue) 
	{
		this.minValue = minValue;
		return this;
	}

	public RangeParameter<T> maxValue(T maxValue) 
	{
		this.maxValue = maxValue;
		return this;
	}
	
	public RangeParameter<T> setMinExclusive()
	{
		this.isMinInclusive = false;
		return this;
	}
	
	public RangeParameter<T> setMaxInclusive()
	{
		this.isMaxInclusive = true;
		return this;
	}
	
	public boolean isValueValid(T value)
	{
		if(value.compareTo(minValue) < 0 || value.compareTo(maxValue) > 0
				|| !isMinInclusive && value.equals(minValue)
				|| !isMaxInclusive && value.equals(maxValue))
			return false;
		else
			return true;
	}
	
	public RangeParameter<T> defaultValue(String value)
	{
		super.setOptional();
		String error = parse(value);
		if(error != null)
		{
			System.err.println("(RangeParameter) Error while parsing the default value: " + error);
			System.exit(-1);
		}
		return this;
	}
	
	public abstract String parse(String value);


	@Override
	public String getValueAsString() {
		return min+","+max;
	}

	public T getMin()
	{
		return min;
	}
	
	public T getMax()
	{
		return max;
	}
}
