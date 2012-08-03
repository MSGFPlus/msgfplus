package edu.ucsd.msjava.msscorer;

import edu.ucsd.msjava.msgf.ScoredSpectrum;
import edu.ucsd.msjava.msutil.Matter;

public interface SimpleDBSearchScorer<T extends Matter> extends ScoredSpectrum<T> {
	// fromIndex: inclusive, toIndex: exclusive
	public int getScore(double[] prefixMassArr, int[] intPrefixMassArr, int fromIndex, int toIndex, int numMods);
}
