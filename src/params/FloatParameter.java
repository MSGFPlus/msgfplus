package params;

public class FloatParameter extends NumberParameter<Float> {

	public FloatParameter(String key, String name, String description) 
	{
		super(key, name, description);
		super.minValue = Float.NEGATIVE_INFINITY;
		super.maxValue = Float.POSITIVE_INFINITY;
	}

	@Override
	public String parse(String value) {
		try {
			super.value = Float.valueOf(value);
			if(this.value < minValue || this.value >= maxValue)
				return "must be in the range [" + minValue + "," + maxValue + ")";
		} catch (NumberFormatException e)
		{
			return "must be a float";
		} 
		return null;
	}
}
