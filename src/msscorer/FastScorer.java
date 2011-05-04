package msscorer;

import java.util.ArrayList;

import msdbsearch.CandidatePeptideTree;
import msdbsearch.DatabaseMatch;
import msgf.NominalMass;
import msgf.NominalMassFactory;

// this does not use edge scores
public class FastScorer implements SimpleDBSearchScorer<NominalMass> {

	protected float[] prefixScore = null;
	protected float[] suffixScore = null;
	
	public FastScorer(NominalMassFactory factory, NewScoredSpectrum<NominalMass> scoredSpec, int peptideMass)
	{
		prefixScore = new float[peptideMass];
		suffixScore = new float[peptideMass];
		for(int i=0; i<prefixScore.length; i++)
			prefixScore[i] = Float.MIN_VALUE;
		for(int nominalMass=1; nominalMass<peptideMass; nominalMass++)
		{
			NominalMass node = factory.getInstanceOfIndex(nominalMass);
			if(node != null)
			{
				prefixScore[nominalMass] = scoredSpec.getNodeScore(node, true);
				suffixScore[nominalMass] = scoredSpec.getNodeScore(node, false);
			}
		}
	}
	
	// fromIndex: inclusive, toIndex: exclusive
	@Override
	public int getScore(double[] prefixMassArr, int[] nominalPrefixMassArr, int fromIndex, int toIndex)
	{
		int score = 0;
		int peptideMass = nominalPrefixMassArr[toIndex-1];
		for(int i=fromIndex; i<toIndex-1; i++)
		{
			int prefixMass = nominalPrefixMassArr[i];
			int suffixMass = peptideMass - prefixMass;
			int curScore = Math.round(prefixScore[prefixMass]+suffixScore[suffixMass]);
			score += curScore;
		}
		return score;
	}

	public ArrayList<DatabaseMatch> getPeptideMatches(CandidatePeptideTree candidates, int scoreThreshold)
	{
		return null;
	}
	
	@Override
	public int getNodeScore(NominalMass prefixMass, NominalMass suffixMass) {
		return Math.round(prefixScore[prefixMass.getNominalMass()]+suffixScore[suffixMass.getNominalMass()]);
	}

	@Override
	public int getEdgeScore(NominalMass curNode, NominalMass prevNode, float theoMass) {
		return 0;
	}
}
