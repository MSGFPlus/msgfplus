package msscorer;

import msgf.ScoredSpectrum;
import msutil.Matter;

public interface SimpleDBSearchScorer<T extends Matter> extends ScoredSpectrum<T> {
	public int getScore(double[] prefixMassArr, int[] intPrefixMassArr, int fromIndex, int toIndex);
}
