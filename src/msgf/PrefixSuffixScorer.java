package msgf;

import msgf.IntMassFactory.IntMass;
import msscorer.NewScoredSpectrum;
import msscorer.SimpleDBSearchScorer;
import msutil.ActivationMethod;
import msutil.Peak;

public class PrefixSuffixScorer implements SimpleDBSearchScorer<IntMass> {

	private float[] prefixScore = null;
	private float[] suffixScore = null;
	private final boolean useVariablePMBin;
	private boolean mainIonDirection;
	private final String activationMethodName;
	private final Peak precursor;
	private final int scanNum;
	
	public PrefixSuffixScorer(IntMassFactory factory, ScoredSpectrum<IntMass> scoredSpec, int peptideMassIndex)
	{
		prefixScore = new float[peptideMassIndex];
		suffixScore = new float[peptideMassIndex];
		for(int massIndex=1; massIndex<peptideMassIndex; massIndex++)
		{
			IntMass node = factory.getInstanceOfIndex(massIndex);
			if(node != null)
			{
				prefixScore[massIndex] = scoredSpec.getNodeScore(node, true);
				suffixScore[massIndex] = scoredSpec.getNodeScore(node, false);
			}
		}
		mainIonDirection = scoredSpec.getMainIonDirection();
		if(Math.round(factory.getRescalingConstant()) == 1)
			useVariablePMBin = true;
		else
			useVariablePMBin = false;
		this.activationMethodName = scoredSpec.getActivationMethodName();
		this.precursor = scoredSpec.getPrecursorPeak();
		this.scanNum = scoredSpec.getScanNum();
	}
	
	// fromIndex: inclusive, toIndex: exclusive
	@Override
	public int getScore(double[] prefixMassArr, int[] intPrefixMassArr, int fromIndex, int toIndex, int numMods)
	{
		int score = 0;
		int theoPeptideMass = intPrefixMassArr[toIndex-1];
		int peptideMass = useVariablePMBin ? theoPeptideMass : prefixScore.length;
		for(int i=fromIndex; i<toIndex-1; i++)
		{
			int suffixMass = theoPeptideMass - intPrefixMassArr[i];
			int prefixMass = peptideMass - suffixMass;
			score += Math.round(prefixScore[prefixMass]+suffixScore[suffixMass]);
		}
		score += FlexAminoAcidGraph.MODIFIED_EDGE_PENALTY*numMods;

		return score;
	}

	@Override
	public int getNodeScore(IntMass prefixMass, IntMass suffixMass) {
		return Math.round(prefixScore[prefixMass.getNominalMass()]+suffixScore[suffixMass.getNominalMass()]);
	}

	@Override
	public int getEdgeScore(IntMass curNode, IntMass prevNode, float edgeMass) {
		assert(false);	// do not use
		return 0;
	}

	@Override
	public boolean getMainIonDirection() {
		return mainIonDirection;
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
	public float getNodeScore(IntMass node, boolean isPrefix) {
		if(isPrefix)
			return prefixScore[node.getNominalMass()];
		else
			return suffixScore[node.getNominalMass()];
	}

	@Override
	public int getScanNum() {
		return scanNum;
	}
}
