package msgf;

import java.util.List;

import msutil.Matter;
import msutil.Peak;

public class ScoredSpectrumSum<T extends Matter> implements ScoredSpectrum<T> {

	private List<ScoredSpectrum<T>> scoredSpecList;
	private final String activationMethodName;
	private Peak precursor;
	private int scanNum;
	
	public ScoredSpectrumSum(List<ScoredSpectrum<T>> scoredSpecList)
	{
		this.scoredSpecList = scoredSpecList;
		StringBuffer buf = null;
		for(ScoredSpectrum<T> scoredSpec : scoredSpecList)
		{
			if(buf != null)
				buf.append("/"+scoredSpec.getActivationMethodName());
			else
			{
				buf = new StringBuffer();
				buf.append(scoredSpec.getActivationMethodName());
				precursor = scoredSpec.getPrecursorPeak().clone();
				scanNum = scoredSpec.getScanNum(); 
			}
		}
		activationMethodName = buf.toString();
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
	public String getActivationMethodName() {
		return activationMethodName;
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
	public int getScanNum() {
		return scanNum;
	}
}
