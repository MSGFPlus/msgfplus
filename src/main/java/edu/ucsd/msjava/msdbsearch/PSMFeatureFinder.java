package edu.ucsd.msjava.msdbsearch;

import java.util.ArrayList;
import java.util.List;

import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.NewScoredSpectrum;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.Pair;

public class PSMFeatureFinder {
	
	private final Spectrum spec;	// MS/MS spectrum
//	private final Spectrum precursorSpec;
	private final Peptide peptide;
	private final NewScoredSpectrum<NominalMass> scoredSpec;
	
	private Float ms2IonCurrent = null;	// summed intensity of all observed product ions
	private Float nTermIonCurrent = null;	// summed intensity of all explained N-term product ions
	private Float cTermIonCurrent = null;	// summed intensity of all explained C-term product ions
	
	private Float errSDAll = null;
	private Float errMeanAll = null;
	private Float errSD7 = null;
	private Float errMean7 = null;
	
	private Float errRSDAll = null;
	private Float errRMeanAll = null;
	private Float errRSD7 = null;
	private Float errRMean7 = null;
	
//	private Float ms1IonCurrent;
//	private Float isolationWindowEfficiency;
	private Tolerance mme;
	
	public PSMFeatureFinder(Spectrum spec, Spectrum precursorSpec, Peptide peptide, NewRankScorer scorer)
	{
		this.spec = spec;
		this.peptide = peptide;
//		this.precursorSpec = precursorSpec;

		scoredSpec = scorer.getScoredSpectrum(spec);
		if(scorer.getSpecDataType().getInstrumentType().isHighResolution())
			mme = new Tolerance(20f, true);	// for high-precision MS/MS, set tolerance as 20ppm
		else
			mme = new Tolerance(0.5f, false);	// low resolution: 0.5Da

		extractFeatures();
	}

	public PSMFeatureFinder(Spectrum spec, Peptide peptide, NewRankScorer scorer)
	{
		this(spec, null, peptide, scorer);
	}
	
	public List<Pair<String, String>> getAllFeatures()
	{
		List<Pair<String, String>> list = new ArrayList<Pair<String, String>>();
		
		Float explainedIonCurrentRatio = getExplainedIonCurrent();
		if(explainedIonCurrentRatio != null)
			list.add(new Pair<String,String>("ExplainedIonCurrentRatio", String.valueOf(getExplainedIonCurrent())));

		Float nTermExplainedIonCurrent = getNTermExplainedIonCurrent();
		if(nTermExplainedIonCurrent != null)
			list.add(new Pair<String,String>("NTermIonCurrentRatio", String.valueOf(nTermExplainedIonCurrent)));

		Float cTermExplainedIonCurrent = getCTermExplainedIonCurrent();
		if(cTermExplainedIonCurrent != null)
			list.add(new Pair<String,String>("CTermIonCurrentRatio", String.valueOf(cTermExplainedIonCurrent)));

		Float ms2IonCurrent = getMS2IonCurrent();
		if(explainedIonCurrentRatio != null)
			list.add(new Pair<String,String>("MS2IonCurrent", String.valueOf(ms2IonCurrent)));

		Float ms1IonCurrent = getMS1IonCurrent();
		if(ms1IonCurrent != null)
			list.add(new Pair<String,String>("MS1IonCurrent", String.valueOf(ms1IonCurrent)));
		
		Float isolationWindowEfficiency = getIsolationWindowEfficiency();
		if(isolationWindowEfficiency != null)
			list.add(new Pair<String,String>("IsolationWindowEfficiency", String.valueOf(isolationWindowEfficiency)));

		if(this.errMeanAll != null)
			list.add(new Pair<String,String>("MeanErrorAll", String.valueOf(errMeanAll)));

		if(this.errSDAll != null)
			list.add(new Pair<String,String>("StdevErrorAll", String.valueOf(errSDAll)));

		if(this.errMean7 != null)
			list.add(new Pair<String,String>("MeanErrorTop7", String.valueOf(errMean7)));

		if(this.errSD7 != null)
			list.add(new Pair<String,String>("StdevErrorTop7", String.valueOf(errSD7)));

		if(this.errRMeanAll != null)
			list.add(new Pair<String,String>("MeanRelErrorAll", String.valueOf(errRMeanAll)));

		if(this.errRSDAll != null)
			list.add(new Pair<String,String>("StdevRelErrorAll", String.valueOf(errRSDAll)));

		if(this.errRMean7 != null)
			list.add(new Pair<String,String>("MeanRelErrorTop7", String.valueOf(errRMean7)));

		if(this.errRSD7 != null)
			list.add(new Pair<String,String>("StdevRelErrorTop7", String.valueOf(errRSD7)));
		
		return list;
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
		
		MassErrorStat errStat = new MassErrorStat();
		double prm = 0, srm = 0;
		for(int i=0; i<peptide.size()-1; i++)
		{
			prm += peptide.get(i).getAccurateMass();
			srm += peptide.get(peptide.size()-1-i).getAccurateMass();
			nTermIonCurrent += scoredSpec.getExplainedIonCurrent((float)prm, true, mme);
			cTermIonCurrent += scoredSpec.getExplainedIonCurrent((float)srm, false, mme);
			
			Pair<Float, Float> err;
			if((err = scoredSpec.getMassErrorWithIntensity((float)prm, true, mme)) != null)
				errStat.add(err);
			if((err = scoredSpec.getMassErrorWithIntensity((float)srm, false, mme)) != null)
				errStat.add(err);
		}

		if(errStat.size() > 0)
		{
			errStat.computeStats();
			this.errMeanAll = errStat.getMean();
			this.errSDAll = errStat.getSd();
			this.errMean7 = errStat.getMean7();
			this.errSD7 = errStat.getSd7();
			
			this.errRMeanAll = errStat.getRMean();
			this.errRSDAll = errStat.getRSd();
			this.errRMean7 = errStat.getRMean7();
			this.errRSD7 = errStat.getRSd7();
		}
		
		this.nTermIonCurrent = nTermIonCurrent;
		this.cTermIonCurrent = cTermIonCurrent;
	}
	
	public Float getExplainedIonCurrent()
	{
		Float nEIC = getNTermExplainedIonCurrent();
		Float cEIC = getCTermExplainedIonCurrent();
		
		if(nEIC != null && cEIC != null)
			return nEIC+cEIC;
		else
			return null;
	}	
	
	public Float getNTermExplainedIonCurrent()
	{
		return nTermIonCurrent/ms2IonCurrent;
	}
	
	public Float getCTermExplainedIonCurrent()
	{
		return cTermIonCurrent/ms2IonCurrent;
	}
	
	public Float getMS2IonCurrent()
	{
		return ms2IonCurrent;
	}
	
	public Float getMS1IonCurrent()
	{
//		return ms1IonCurrent;
		return null;
	}
	
	public Float getIsolationWindowEfficiency()
	{
//		return isolationWindowEfficiency;
		return null;
	}
}
