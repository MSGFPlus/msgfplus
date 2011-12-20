package params;

import msgf.Tolerance;

public class ToleranceParameter extends Parameter {

	private Tolerance leftTolerance;
	private Tolerance rightTolerance;
	
	public ToleranceParameter(String key, String name, String description) {
		super(key, name, description);
	}

	public ToleranceParameter defaultValue(String value)
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
		if(token.length == 1)
		{
			leftTolerance = rightTolerance = Tolerance.parseToleranceStr(token[0]);
		}
		else if(token.length == 2)
		{
			leftTolerance = Tolerance.parseToleranceStr(token[0]);
			rightTolerance = Tolerance.parseToleranceStr(token[1]);
		}
		if(leftTolerance == null || rightTolerance == null)
		{
			return "illegal tolerance value";
		}
		if(leftTolerance.isTolerancePPM() != rightTolerance.isTolerancePPM())
		{
			return "left and right tolerance units must be the same";
		}
		if(leftTolerance.getValue() < 0 || rightTolerance.getValue() < 0)
		{
			return "parent mass tolerance must not be negative";
		}
		return null;
	}

	@Override
	public String getValueAsString() {
		if(leftTolerance == null || rightTolerance == null)
			return null;
		return leftTolerance.toString()+","+rightTolerance.toString();
	}

	public Tolerance getLeftTolerance()
	{
		return leftTolerance;
	}
	
	public Tolerance getRightTolerance()
	{
		return rightTolerance;
	}
}
