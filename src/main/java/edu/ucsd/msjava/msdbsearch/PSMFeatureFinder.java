package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.NewScoredSpectrum;
import edu.ucsd.msjava.msscorer.NewScorerFactory;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msutil.IonType;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.Spectrum;

public class PSMFeatureFinder {
	
	private final Spectrum spec;	// MS/MS spectrum
	private final Spectrum precursorSpec;
	private final Peptide peptide;
	private final NewScoredSpectrum<NominalMass> scoredSpec;
	
	private Float ms2IonCurrent = null;	// summed intensity of all observed product ions
	private Float nTermIonCurrent = null;	// summed intensity of all explained N-term product ions
	private Float cTermIonCurrent = null;	// summed intensity of all explained C-term product ions
	
	private Float ms1IonCurrent;
	private Float isolationWindowEfficiency;
	
	public PSMFeatureFinder(Spectrum spec, Peptide peptide, Spectrum precursorSpec, SpecDataType specDataType)
	{
		this.spec = spec;
		this.peptide = peptide;
		this.precursorSpec = precursorSpec;
		NewRankScorer scorer = NewScorerFactory.get(spec.getActivationMethod(), specDataType.getInstrumentType(), specDataType.getEnzyme(), specDataType.getProtocol());
		scoredSpec = scorer.getScoredSpectrum(spec);
		
		extractFeatures();
	}

	public PSMFeatureFinder(Spectrum spec, Peptide peptide, SpecDataType specDataType)
	{
		this(spec, peptide, null, specDataType);
	}
	
	private void extractFeatures()
	{
		computeSumIonCurrent();
		computeExplainedIonCurrent();
	}
	
	private void computeSumIonCurrent()
	{
		float ms2IonCurrent = 0f;
		for(Peak p : spec)
			ms2IonCurrent += p.getIntensity();
		
		this.ms2IonCurrent = ms2IonCurrent;
	}
	
	private void computeExplainedIonCurrent()
	{
		float nTermIonCurrent = 0f, cTermIonCurrent = 0f;
		
		double prm = 0, srm = 0;
		for(int i=0; i<peptide.size()-1; i++)
		{
			prm += peptide.get(i).getAccurateMass();
			srm += peptide.get(peptide.size()-1-i).getAccurateMass();
			nTermIonCurrent += scoredSpec.getExplainedIonCurrent((float)prm, true);
			cTermIonCurrent += scoredSpec.getExplainedIonCurrent((float)srm, false);
		}
		
		this.nTermIonCurrent = nTermIonCurrent;
		this.cTermIonCurrent = cTermIonCurrent;
	}
	
	public Float getExplainedIonCurrent()
	{
		Float nEIC = getNTermExplainedIonCurrent();
		Float cEIC = getNTermExplainedIonCurrent();
		
		if(nEIC != null && cEIC != null)
			return nEIC+cEIC;
		else
			return null;
	}	
	
	public Float getNTermExplainedIonCurrent()
	{
		return nTermIonCurrent/ms2IonCurrent;
	}
	
	public float getCTermExplainedIonCurrent()
	{
		return cTermIonCurrent/ms2IonCurrent;
	}
	
	public float getMS1IonCurrent()
	{
		return ms1IonCurrent;
	}
	
	public float getMS2IonCurrent()
	{
		return ms2IonCurrent;
	}
	
	public float getIsolationWindowEfficiency()
	{
		return isolationWindowEfficiency;
	}
}
