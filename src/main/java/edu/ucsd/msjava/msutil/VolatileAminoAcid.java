package edu.ucsd.msjava.msutil;

import java.util.HashMap;

public class VolatileAminoAcid extends AminoAcid {

	private VolatileAminoAcid(float mass) {
		super('*', String.format("(%.3f)", mass), mass);
	}
	
	@Override
	public String getResidueStr()
	{
		return super.getName();
	}	

	@Override
	public boolean isModified()
	{
		return true;
	}
	
	public static AminoAcid getVolatileAminoAcid(float mass)
	{
		AminoAcid aa = table.get(mass);
		if(aa == null)
		{
//			System.out.println("Register " + mass);
			aa = new VolatileAminoAcid(mass);
			table.put(mass, aa);
		}
		return aa;
	}
	
	private static HashMap<Float,AminoAcid> table = new HashMap<Float,AminoAcid>();
}
