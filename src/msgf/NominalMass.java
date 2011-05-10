package msgf;

import msutil.Constants;
import msutil.Matter;

public class NominalMass extends Matter {
	private int nominalMass;
	public NominalMass(int nominalMass)
	{
		this.nominalMass = nominalMass;
	}
	@Override
	public float getMass() {
		return nominalMass/Constants.INTEGER_MASS_SCALER;
	}

	@Override
	public int getNominalMass() {
		return nominalMass;
	}

	@Override
	public int hashCode()	
	{ 
		return nominalMass; 
	}
	
	@Override
	public boolean equals(Object obj)	
	{
		if(!(obj instanceof NominalMass))
			return false;
		return (nominalMass == ((NominalMass)obj).nominalMass);
	}
	
	@Override
	public String toString()
	{
		return String.valueOf(nominalMass);
	}
	
	public static int toNominalMass(float mass)
	{
		return Math.round(mass*Constants.INTEGER_MASS_SCALER);
	}
	
	public static float getMassFromNominalMass(int nominalMass)
	{
		return nominalMass/Constants.INTEGER_MASS_SCALER;
	}
}