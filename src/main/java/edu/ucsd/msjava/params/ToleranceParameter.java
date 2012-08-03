package edu.ucsd.msjava.params;

import edu.ucsd.msjava.msgf.Tolerance;

public class ToleranceParameter extends Parameter {

	private Tolerance leftTolerance;
	private Tolerance rightTolerance;
	private boolean allowAsymmetricValues = true;
	
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
	
	public ToleranceParameter doNotAllowAsymmetricValues()
	{
		this.allowAsymmetricValues = false;
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
			if(allowAsymmetricValues)
			{
				leftTolerance = Tolerance.parseToleranceStr(token[0]);
				rightTolerance = Tolerance.parseToleranceStr(token[1]);
			}
			else
				return "asymmetric values are not allowed";
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
