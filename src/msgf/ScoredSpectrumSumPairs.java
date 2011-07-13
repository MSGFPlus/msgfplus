package msgf;

import msutil.Matter;
import msutil.Peak;

public class ScoredSpectrumSumPairs<T extends Matter> implements ScoredSpectrum<T> {

	private ScoredSpectrum<T> scoredSpec1;
	private ScoredSpectrum<T> scoredSpec2;
	private String activationMethodName;
	
	public ScoredSpectrumSumPairs(ScoredSpectrum<T> scoredSpec1, ScoredSpectrum<T> scoredSpec2)
	{
		this.scoredSpec1 = scoredSpec1;
		this.scoredSpec2 = scoredSpec2;
		activationMethodName = scoredSpec1.getActivationMethodName()+"/"+scoredSpec2.getActivationMethodName();
	}
	
	@Override
	public int getNodeScore(T prefixResidueNode, T suffixResidueNode) {
		return scoredSpec1.getNodeScore(prefixResidueNode, suffixResidueNode)+scoredSpec2.getNodeScore(prefixResidueNode, suffixResidueNode);
	}

	@Override
	public int getEdgeScore(T curNode, T prevNode, float theoMass) {
		return scoredSpec1.getEdgeScore(curNode, prevNode, theoMass)+scoredSpec2.getEdgeScore(curNode, prevNode, theoMass);
	}
	
	@Override
	public boolean getMainIonDirection() {
		assert(false): "Not supported!";
		return false;
	}

	@Override
	public String getActivationMethodName() {
		return activationMethodName;
	}

	@Override
	public Peak getPrecursorPeak() {
		return scoredSpec1.getPrecursorPeak();
	}

	@Override
	public float getNodeScore(T node, boolean isPrefix) {
		return scoredSpec1.getNodeScore(node, isPrefix)+scoredSpec2.getNodeScore(node, isPrefix);
	}
	
	@Override
	public int getScanNum() {
		return scoredSpec1.getScanNum();
	}
	
}
