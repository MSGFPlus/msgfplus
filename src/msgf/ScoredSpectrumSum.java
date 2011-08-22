package msgf;

import java.util.List;

import msutil.ActivationMethod;
import msutil.Matter;
import msutil.Peak;

public class ScoredSpectrumSum<T extends Matter> implements ScoredSpectrum<T> {

	private List<ScoredSpectrum<T>> scoredSpecList;
	private Peak precursor;
	private ActivationMethod[] activationMethodArr;
	private int[] scanNumArr;
	
	public ScoredSpectrumSum(List<ScoredSpectrum<T>> scoredSpecList)
	{
		this.scoredSpecList = scoredSpecList;
		scanNumArr = new int[scoredSpecList.size()];
		activationMethodArr = new ActivationMethod[scoredSpecList.size()];

		int i=0;
		for(ScoredSpectrum<T> scoredSpec : scoredSpecList)
		{
			scanNumArr[i] = scoredSpec.getScanNumArr()[0];
			activationMethodArr[i] = scoredSpec.getActivationMethodArr()[0];
		}
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

	@Override
	public Peak getPrecursorPeak() {
		return precursor;
	}

	@Override
	public float getNodeScore(T node, boolean isPrefix) {
		float sum = 0;
		for(ScoredSpectrum<T> scoredSpec : scoredSpecList)
			sum += scoredSpec.getNodeScore(node, isPrefix);
		return sum;
	}

	@Override
	public ActivationMethod[] getActivationMethodArr() {
		return this.activationMethodArr;
	}

	@Override
	public int[] getScanNumArr() {
		return this.scanNumArr;
	}
}
