package edu.ucsd.msjava.msutil;

public class ScanType
{
	public ScanType(ActivationMethod activationMethod, boolean isHighPrecision, int msLevel) 
	{
		this.activationMethod = activationMethod;
		this.msLevel = msLevel;
		this.isHighPrecision = isHighPrecision;
	}
    
	public ScanType(ActivationMethod activationMethod, boolean isHighPrecision, int msLevel, float scanStartTime) 
	{
		this.activationMethod = activationMethod;
		this.msLevel = msLevel;
		this.isHighPrecision = isHighPrecision;
        this.scanStartTime = scanStartTime;
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
    public float getScanStartTime()
    {
        return scanStartTime;
    }
	
	private ActivationMethod activationMethod;
	private int msLevel;
	private boolean isHighPrecision;
    private float scanStartTime;
}

