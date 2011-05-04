package msgf;

import java.io.Serializable;

public class Tolerance implements Serializable { // Serializable is needed in order to make RankScorer serializable
	/**
	 * 
	 */
	public static final Tolerance ZERO_TOLERANCE = new Tolerance(0);
	
	private static final long serialVersionUID = 1L;
	private float value;
	private boolean isTolerancePPM = false;
	
	public Tolerance(float value)
	{
		this(value, false);
	}
	
	public Tolerance(float value, boolean isTolerancePPM)
	{
		this.value = value;
		this.isTolerancePPM = isTolerancePPM;
	}
	
	public float getValue()			{ return value; }
	public boolean isTolerancePPM()	{ return isTolerancePPM; }
	public float getToleranceAsDa(float mass)
	{
		if(!isTolerancePPM)
			return value;
		else
			return 1e-6f*value*mass;
	}
	
	// added by Kyowon
	public float getToleranceAsPPM(float mass)
	{
		if(isTolerancePPM)
			return value;
		else return value * 1e6f / mass;
	}
	
	public String toString()
	{
		if(!isTolerancePPM)
			return value+"Da";
		else
			return value+"ppm";
	}
	
	public static Tolerance parseToleranceStr(String tolStr)
	{
		Float val = null;
		boolean isTolerancePPM = false;
		if(tolStr.endsWith("Da"))
		{
			try {
				val = Float.parseFloat(tolStr.substring(0, tolStr.length()-2).trim()); 
				isTolerancePPM = false;
			}
			catch (NumberFormatException e) {}
		}
		else if(tolStr.endsWith("ppm"))
		{
			try {
				val = Float.parseFloat(tolStr.substring(0, tolStr.length()-3).trim()); 
				isTolerancePPM = true;
			}
			catch (NumberFormatException e) {}
		}
		if(val == null)
			return null;
		else
			return new Tolerance(val, isTolerancePPM);
		
	}
}
