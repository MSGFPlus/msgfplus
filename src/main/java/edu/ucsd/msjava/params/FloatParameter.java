package edu.ucsd.msjava.params;

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
			String range = (super.isMinInclusive ? "[" : "(") + minValue + "," + maxValue + (super.isMaxInclusive ? "]" : ")"); 
			if(this.value < minValue || this.value > maxValue
				|| !super.isMinInclusive && this.value.equals(minValue)
				|| !super.isMaxInclusive && this.value.equals(maxValue))
				return "must be in the range " + range;
		} catch (NumberFormatException e)
		{
			return "must be a float";
		} 
		return null;
	}
}
