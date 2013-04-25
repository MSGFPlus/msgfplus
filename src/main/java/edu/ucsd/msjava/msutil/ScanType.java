package edu.ucsd.msjava.msutil;

public class ScanType
{
	public ScanType(ActivationMethod activationMethod, boolean isHighPrecision, int msLevel) 
	{
		this.activationMethod = activationMethod;
		this.msLevel = msLevel;
		this.isHighPrecision = isHighPrecision;
	}
	
	public ActivationMethod getActivationMethod() 
	{
		return activationMethod;
	}
	public int getMsLevel() 
	{
		return msLevel;
	}
	public boolean isHighPrecision() 
	{
		return isHighPrecision;
	}
	
	private ActivationMethod activationMethod;
	private int msLevel;
	private boolean isHighPrecision;
}

