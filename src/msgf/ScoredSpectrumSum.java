package msgf;

import java.util.List;

import msutil.Matter;

public class ScoredSpectrumSum<T extends Matter> implements ScoredSpectrum<T> {

	private List<ScoredSpectrum<T>> scoredSpecList;
	
	public ScoredSpectrumSum(List<ScoredSpectrum<T>> scoredSpecList)
	{
		this.scoredSpecList = scoredSpecList;
	}
	
	@Override
	public int getNodeScore(T prefixResidueNode, T suffixResidueNode) {
		int sum = 0;
		for(ScoredSpectrum<T> scoredSpec : scoredSpecList)
			sum += scoredSpec.getNodeScore(prefixResidueNode, suffixResidueNode);
		return sum;
	}

	@Override
	public int getEdgeScore(T curNode, T prevNode, float theoMass) {
		int sum = 0;
		for(ScoredSpectrum<T> scoredSpec : scoredSpecList)
			sum += scoredSpec.getEdgeScore(curNode, prevNode, theoMass);
		return sum;
	}

	@Override
	public boolean getMainIonDirection() {
		assert(false): "Not supported!";
		return false;
	}
}
