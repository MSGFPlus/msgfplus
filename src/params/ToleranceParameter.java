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
		if(parse(value) == false)
		{
			System.err.println("(ToleranceParameter) Error while setting default value: " + value);
			System.exit(-1);
		}
		return this;
	}
	
	@Override
	public boolean parse(String value) {
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
//			printUsageAndExit("Illegal tolerance value: " + value);
			return false;
		}
		if(leftTolerance.isTolerancePPM() != rightTolerance.isTolerancePPM())
		{
//			printUsageAndExit("Left and right tolerance units must be the same: " + value);
			return false;
		}
		if(leftTolerance.getValue() < 0 || rightTolerance.getValue() < 0)
		{
//			printUsageAndExit("Parent mass tolerance must not be negative: " + value);
			return false;
		}
		return true;
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
