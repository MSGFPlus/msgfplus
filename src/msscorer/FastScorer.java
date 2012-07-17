package msscorer;

import msgf.FlexAminoAcidGraph;
import msgf.NominalMass;
import msgf.ScoredSpectrum;
import msutil.ActivationMethod;
import msutil.Composition;
import msutil.Peak;

// this does not use edge scores
public class FastScorer implements SimpleDBSearchScorer<NominalMass> {

	protected float[] prefixScore = null;
	protected float[] suffixScore = null;
	private boolean mainIonDirection;
	protected Peak precursor;
	protected ActivationMethod[] activationMethodArr;
	private int[] scanNumArr;
	
	public FastScorer(ScoredSpectrum<NominalMass> scoredSpec, int peptideMass)
	{
		prefixScore = new float[peptideMass];
		suffixScore = new float[peptideMass];
		for(int i=0; i<prefixScore.length; i++)
			prefixScore[i] = Float.MIN_VALUE;
		for(int nominalMass=1; nominalMass<peptideMass; nominalMass++)
		{
			NominalMass node = new NominalMass(nominalMass);
			prefixScore[nominalMass] = scoredSpec.getNodeScore(node, true);
			suffixScore[nominalMass] = scoredSpec.getNodeScore(node, false);
		}
		mainIonDirection = scoredSpec.getMainIonDirection();
		
		this.precursor = scoredSpec.getPrecursorPeak();
		this.activationMethodArr = scoredSpec.getActivationMethodArr();
		this.scanNumArr = scoredSpec.getScanNumArr();
	}

	public Peak getPrecursorPeak()	{ return precursor; }
	
	public ActivationMethod[] getActivationMethodArr()	{ return activationMethodArr; }
	
	public float getParentMass()	{ return precursor.getMass(); }
	public float getPeptideMass()	{ return precursor.getMass()-(float)(Composition.H2O); }
	public int getCharge()			{ return precursor.getCharge(); }
	
	
	// fromIndex: inclusive, toIndex: exclusive
	public int getScore(double[] prefixMassArr, int[] nominalPrefixMassArr, int fromIndex, int toIndex, int numMods)
	{
		int score = 0;
		int peptideMass = nominalPrefixMassArr[toIndex-1];
		for(int i=fromIndex; i<toIndex-1; i++)
		{
			int prefixMass = nominalPrefixMassArr[i];
			int suffixMass = peptideMass - prefixMass;
//			if(prefixMass >= prefixScore.length || suffixMass >= suffixScore.length)
//			{
//				System.out.println("Debug");
//			}
			int curScore = Math.round(prefixScore[prefixMass]+suffixScore[suffixMass]);
			score += curScore;
		}
		
		score += FlexAminoAcidGraph.MODIFIED_EDGE_PENALTY*numMods;
		return score;
	}

	public int getNodeScore(NominalMass prefixMass, NominalMass suffixMass) {
		if(prefixMass.getNominalMass() >= prefixScore.length ||
				suffixMass.getNominalMass() >= suffixScore.length)
			System.out.println("Debug");
		return Math.round(prefixScore[prefixMass.getNominalMass()]+suffixScore[suffixMass.getNominalMass()]);
	}

	public int getEdgeScore(NominalMass curNode, NominalMass prevNode, float theoMass) {
		return 0;
	}
	
	public boolean getMainIonDirection() {
		return mainIonDirection;
	}

	public float getNodeScore(NominalMass node, boolean isPrefix) {
		if(isPrefix)
			return prefixScore[node.getNominalMass()];
		else
			return suffixScore[node.getNominalMass()];
	}

	public int[] getScanNumArr() {
		return scanNumArr;
	}
	
}
