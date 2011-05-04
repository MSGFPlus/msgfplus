package msgf;

import msgf.IntMassFactory.IntMass;
import msscorer.NewScoredSpectrum;
import msscorer.SimpleDBSearchScorer;

public class PrefixSuffixScorer implements SimpleDBSearchScorer<IntMass> {

	private float[] prefixScore = null;
	private float[] suffixScore = null;
	private final boolean useVariablePMBin;
	
	public PrefixSuffixScorer(IntMassFactory factory, NewScoredSpectrum<IntMass> scoredSpec, int peptideMassIndex)
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
		if(Math.round(factory.getRescalingConstant()) == 1)
			useVariablePMBin = true;
		else
			useVariablePMBin = false;
	}
	
	// fromIndex: inclusive, toIndex: exclusive
	public int getScore(double[] prefixMassArr, int[] intPrefixMassArr, int fromIndex, int toIndex)
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
}
