package edu.ucsd.msjava.params;

public class FloatRangeParameter extends RangeParameter<Float> {
	public FloatRangeParameter(String key, String name, String description) {
		super(key, name, description);
		super.minValue = Float.MIN_VALUE;
		super.maxValue = Float.MAX_VALUE;
		super.isMinInclusive = true;
		super.isMaxInclusive = false;
	}

	@Override
	public String parse(String value) {
		String[] token = value.split(",");
		try {
			if(token.length == 2)
			{
				min = Float.parseFloat(token[0]);
				max = Float.parseFloat(token[1]);
			}
			else
			{
				return "illegar syntax";
			}
		} catch (NumberFormatException e)
		{
			return "not a valid float or float range";
		} 
		
		if(min >= max || !isValueValid(min) || !isValueValid(max))
		{
			return "not a valid range";
		}
		return null;
	}
}
