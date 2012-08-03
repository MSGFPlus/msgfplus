package edu.ucsd.msjava.msgf;

import java.util.List;

import edu.ucsd.msjava.msutil.Matter;

public class ScoredSpectrumSumHP<T extends Matter> extends ScoredSpectrumSum<T> {
	public ScoredSpectrumSumHP(List<ScoredSpectrum<T>> scoredSpecList)
	{
		super(scoredSpecList);
	}
	
	
}
