package msscorer;

import msgf.NominalMass;
import msgf.ScoredSpectrum;
import msutil.Composition;
import msutil.Peak;

// this does not use edge scores
public class FastScorer implements SimpleDBSearchScorer<NominalMass> {

	protected float[] prefixScore = null;
	protected float[] suffixScore = null;
	private boolean mainIonDirection;
	protected Peak precursor;
	protected String activationMethodName;
	private int scanNum;
	
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
		this.activationMethodName = scoredSpec.getActivationMethodName();
		this.scanNum = scoredSpec.getScanNum();
	}

	@Override
	public Peak getPrecursorPeak()	{ return precursor; }
	
	@Override
	public String getActivationMethodName()	{ return activationMethodName; }
	
	public float getParentMass()	{ return precursor.getMass(); }
	public float getPeptideMass()	{ return precursor.getMass()-(float)(Composition.H2O); }
	public int getCharge()			{ return precursor.getCharge(); }
	
	
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
//			try {
			int curScore = Math.round(prefixScore[prefixMass]+suffixScore[suffixMass]);
			score += curScore;
//			} catch (ArrayIndexOutOfBoundsException e)
//			{
////				System.out.println("*******"+prefixMass+" "+suffixMass+" "+fromIndex+" "+toIndex+" "+scanNum);
//				e.printStackTrace();
//				System.exit(-1);
//			}
		}
		return score;
	}

	@Override
	public int getNodeScore(NominalMass prefixMass, NominalMass suffixMass) {
		if(prefixMass.getNominalMass() >= prefixScore.length ||
				suffixMass.getNominalMass() >= suffixScore.length)
			System.out.println("Debug");
		return Math.round(prefixScore[prefixMass.getNominalMass()]+suffixScore[suffixMass.getNominalMass()]);
	}

	@Override
	public int getEdgeScore(NominalMass curNode, NominalMass prevNode, float theoMass) {
		return 0;
	}
	
	@Override
	public boolean getMainIonDirection() {
		return mainIonDirection;
	}

	@Override
	public float getNodeScore(NominalMass node, boolean isPrefix) {
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
