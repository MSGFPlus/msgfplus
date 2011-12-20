package params;

public class IntRangeParameter extends Parameter {
	private int min = Integer.MIN_VALUE;	// inclusive
	private int max = Integer.MAX_VALUE;	// exclusive

	private int minValue = 0;
	private int maxValue = Integer.MAX_VALUE;
	
	public IntRangeParameter(String key, String name, String description) {
		super(key, name, description);
	}

	public IntRangeParameter minValue(int minValue) 
	{
		this.minValue = minValue;
		return this;
	}

	public IntRangeParameter maxValue(int maxValue) 
	{
		this.maxValue = maxValue;
		return this;
	}
	
	public IntRangeParameter defaultValue(String value)
	{
		super.setOptional();
		String error = parse(value);
		if(error != null)
		{
			System.err.println("(ToleranceParameter) Error while setting default value: " + error);
			System.exit(-1);
		}
		return this;
	}
	
	@Override
	public String parse(String value) {
		String[] token = value.split(",");
		try {
			if(token.length == 1)
			{
				min = Integer.parseInt(token[0]);
				max = min+1;
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
		
		if(min >= max || min < minValue || max >= maxValue)
		{
			return "not a valid range";
		}
		return null;
	}

	@Override
	public String getValueAsString() {
		return min+","+max;
	}

	public int getMin()
	{
		return min;
	}
	
	public int getMax()
	{
		return max;
	}
}
