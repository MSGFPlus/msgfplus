package edu.ucsd.msjava.ipa;

public class Feature {
	public Feature(float mz, int charge, int ms2Scan) 
	{
		this.mz = mz;
		this.charge = charge;
		this.ms2Scan = ms2Scan;
	}
	
	float getMz() 
	{
		return mz;
	}
	int getCharge() 
	{
		return charge;
	}
	int getMS2Scan() 
	{
		return ms2Scan;
	}
	
	@Override
	public boolean equals(Object o)
	{
		if(o instanceof Feature)
		{
			Feature other = (Feature)o;
			if(mz == other.mz 
				&& charge == other.charge 
				&& ms2Scan == other.ms2Scan)
				return true;
		}
		return false;
	}
	
	private float mz;
	private int charge;
	private int ms2Scan;
}
