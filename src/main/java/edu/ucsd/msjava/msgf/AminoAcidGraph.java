package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msutil.Constants;

public class AminoAcidGraph extends GenericDeNovoGraph<NominalMass> {
	
	public AminoAcidGraph(NominalMassFactory factory, float parentMass, ScoredSpectrum<NominalMass> scoredSpec) 
	{
		super(factory, parentMass, Tolerance.ZERO_TOLERANCE, factory.getEnzyme(), scoredSpec);
	}

	public AminoAcidGraph(NominalMassFactory factory, int nominalPeptideMass, ScoredSpectrum<NominalMass> scoredSpec) 
	{
		this(factory, (nominalPeptideMass+18)/Constants.INTEGER_MASS_SCALER, scoredSpec);
	}
	
}
