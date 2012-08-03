package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.Constants;
import edu.ucsd.msjava.msutil.Matter;

public class SpectrumGraphNode extends Matter {
	int index;	// index in the spectrum
	float mass;
	
	public SpectrumGraphNode(int index, float mass)
	{
		this.index = index;
		this.mass = mass;
	}
	
	public int hashCode()
	{
		return index;
	}

	@Override
	public boolean equals(Object obj) {
		if(obj instanceof SpectrumGraphNode)
		{
			if(index == ((SpectrumGraphNode)obj).index)
				return true;
		}
		return false;
	}

	@Override
	public float getMass() {
		return mass;
	}

	@Override
	public int getNominalMass() {
		return Math.round(mass*Constants.INTEGER_MASS_SCALER);
	}
}
