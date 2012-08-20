package edu.ucsd.msjava.params;

public class IntRangeParameter extends RangeParameter<Integer> {
	public IntRangeParameter(String key, String name, String description) {
		super(key, name, description);
		super.minValue = Integer.MIN_VALUE;
		super.maxValue = Integer.MAX_VALUE;
		super.isMinInclusive = true;
		super.isMaxInclusive = false;
	}

	@Override
	public String parse(String value) {
		String[] token = value.split(",");
		try {
			if(token.length == 1)
			{
				min = Integer.parseInt(token[0]);
				max = min;
			}
			else if(token.length == 2)
			{
				min = Integer.parseInt(token[0]);
				max = Integer.parseInt(token[1]);
			}
			else
			{
				return "illegar syntax";
			}
		} catch (NumberFormatException e)
		{
			return "not a valid integer or integer range";
		} 
		
		int minNumber = isMinInclusive ? min : min+1;
		int maxNumber = isMaxInclusive ? max : max-1;
		
		if(minNumber > maxNumber)
		{
			return "not a valid range";
		}
		
//		if(value.compareTo(minValue) < 0 || value.compareTo(maxValue) > 0
//				|| !isMinInclusive && value.equals(minValue)
//				|| !isMaxInclusive && value.equals(maxValue))
//			return false;
		
		return null;
	}
}
