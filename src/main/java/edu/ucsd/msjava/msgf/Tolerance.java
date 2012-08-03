package edu.ucsd.msjava.msgf;

import java.io.Serializable;

public class Tolerance implements Serializable { // Serializable is needed in order to make RankScorer serializable
	/**
	 * 
	 */
	public static final Tolerance ZERO_TOLERANCE = new Tolerance(0);
	
	private static final long serialVersionUID = 1L;
	private float value;
	
	public static enum Unit {
		Da,
		Th,
		PPM,
	}
	
	private final Unit unit;
	
	public Tolerance(float value)
	{
		this(value, false);
	}
	
	// This constructor supports only Da and PPM
	public Tolerance(float value, boolean isTolerancePPM)
	{
		this.value = value;
		if(isTolerancePPM == false)
			unit = Unit.Da;
		else
			unit = Unit.PPM;
	}
	
	public Tolerance(float value, Unit unit)
	{
		this.value = value;
		this.unit = unit;
	}
	
	public float getValue()			{ return value; }
	public Unit getUnit()			{ return unit; }
	
	/**
	 * @return
	 */
	public boolean isTolerancePPM()	{ return unit == Unit.PPM; }
	
	/**
	 * @param mass
	 * @return
	 */
	public float getToleranceAsDa(float mass)
	{
		if(unit == Unit.Th)
		{
			System.err.println("Use getToleranceAsDa(float mass, int charge) instead!");
			System.exit(-1);
		}
		return getToleranceAsDa(mass, 0);
	}
	
	public float getToleranceAsDa(float mass, int charge)
	{
		if(unit == Unit.Da)
			return value;
		else if(unit == Unit.Th)
			return value*charge;
		else 
			return 1e-6f*value*mass;
	}
	
	// added by Kyowon
	public float getToleranceAsPPM(float mass)
	{
		if(unit == Unit.Da)
			return value;
		else return value * 1e6f / mass;
	}
	
	public String toString()
	{
		if(unit == Unit.Da)
			return value+"Da";
		else if(unit == Unit.PPM)
			return value+"ppm";
		else if(unit == Unit.Th)
			return value+"Th";
		else
			return null;
	}
	
	public static Tolerance parseToleranceStr(String tolStr)
	{
		Float val = null;
		Unit unit = null;
		if(tolStr.endsWith("ppm"))
		{
			try {
				val = Float.parseFloat(tolStr.substring(0, tolStr.length()-3).trim()); 
				unit = Unit.PPM;
			}
			catch (NumberFormatException e) {}
		}
		else if(tolStr.endsWith("Da"))
		{
			try {
				val = Float.parseFloat(tolStr.substring(0, tolStr.length()-2).trim()); 
				unit = Unit.Da;
			}
			catch (NumberFormatException e) {}
		}
		else if(tolStr.endsWith("Th"))
		{
			try {
				val = Float.parseFloat(tolStr.substring(0, tolStr.length()-2).trim()); 
				unit = Unit.Th;
			}
			catch (NumberFormatException e) {}
		}
		else
		{
			try {
				val = Float.parseFloat(tolStr); 
				unit = Unit.Da;
			}
			catch (NumberFormatException e) {}
		}
		if(val == null)
			return null;
		else
			return new Tolerance(val, unit);
	}
}
