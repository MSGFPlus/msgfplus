package msgf;

import msutil.Matter;

public class ScoredSpectrumSumPairs<T extends Matter> implements ScoredSpectrum<T> {

	private ScoredSpectrum<T> scoredSpec1;
	private ScoredSpectrum<T> scoredSpec2;
	
	public ScoredSpectrumSumPairs(ScoredSpectrum<T> scoredSpec1, ScoredSpectrum<T> scoredSpec2)
	{
		this.scoredSpec1 = scoredSpec1;
		this.scoredSpec2 = scoredSpec2;
	}
	
	@Override
	public int getNodeScore(T prefixResidueNode, T suffixResidueNode) {
		return scoredSpec1.getNodeScore(prefixResidueNode, suffixResidueNode)+scoredSpec2.getNodeScore(prefixResidueNode, suffixResidueNode);
	}

	@Override
	public int getEdgeScore(T curNode, T prevNode, float theoMass) {
		return scoredSpec1.getEdgeScore(curNode, prevNode, theoMass)+scoredSpec2.getEdgeScore(curNode, prevNode, theoMass);
	}
}
