package msscorer;

import msgf.ScoredSpectrum;
import msutil.Matter;

public interface SimpleDBSearchScorer<T extends Matter> extends ScoredSpectrum<T> {
	// fromIndex: inclusive, toIndex: exclusive
	public int getScore(double[] prefixMassArr, int[] intPrefixMassArr, int fromIndex, int toIndex, int numMods);
}
