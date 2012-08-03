package edu.ucsd.msjava.mslibsearch;

import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Spectrum;

public class ProcessedSpectrum {
	private final Spectrum expSpec;
	private final Spectrum libSpec;
	
	public ProcessedSpectrum(Spectrum expSpec, Spectrum libSpec)
	{
		this.expSpec = expSpec;
		this.libSpec = libSpec;
	}
	
	public Spectrum getSpectrum()
	{
//		boolean[] expPeak = new boolean[NominalMass.toNominalMass(expSpec.getParentMass())];
//		for(Peak p : libSpec)
//		{
//			int nominalMass = NominalMass.toNominalMass(p.getMz());
//			if(nominalMass >= 0 && nominalMass < expPeak.length)
//				expPeak[nominalMass] = true;
//		}
//		
//		Spectrum spec = expSpec.getCloneWithoutPeakList();
//		for(Peak p : expSpec)
//		{
//			int nominalMass = NominalMass.toNominalMass(p.getMz());
//			if(nominalMass >= 0 && nominalMass < expPeak.length && expPeak[NominalMass.toNominalMass(p.getMz())])
//				spec.add(p);
//		}
//		return spec;
		return expSpec;
	}
}
