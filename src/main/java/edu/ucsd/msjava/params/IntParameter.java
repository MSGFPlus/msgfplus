package edu.ucsd.msjava.params;

public class IntParameter extends NumberParameter<Integer> {

	public IntParameter(String key, String name, String description) {
		super(key, name, description);
		super.minValue = 0;
		super.maxValue = Integer.MAX_VALUE;
	}
	
	@Override
	public String parse(String value) {
		try {
			super.value = Integer.valueOf(value);
			String range = (super.isMinInclusive ? "[" : "(") + minValue + "," + maxValue + (super.isMaxInclusive ? "]" : ")"); 
			if(this.value < minValue || this.value > maxValue
					|| !super.isMinInclusive && this.value.equals(minValue)
					|| !super.isMaxInclusive && this.value.equals(maxValue))
				return "must be in the range " + range;
		} catch (NumberFormatException e)
		{
			return "must be an integer";
		} 
		return null;
	}
}
