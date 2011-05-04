package msgf;

public class ScoreDist extends ScoreBound {
	private float[] numDistribution;
	private float[] probDistribution;
	
    ScoreDist(int minScore, int maxScore, boolean calcNumber, boolean calcProb)
	{
		super(minScore, maxScore);
		if(calcNumber)
			numDistribution = new float[maxScore - minScore];
		if(calcProb)
			probDistribution = new float[maxScore - minScore];	
	}
	
	public boolean isProbSet()	{ return probDistribution != null; }
	public boolean isNumSet()	{ return numDistribution != null; }
	
	public void setNumber(int score, float number)	{ numDistribution[score-minScore] = number; }
	public void setProb(int score, float prob)		{ probDistribution[score-minScore] = prob; }

	public void addNumber(int score, float number)	{ numDistribution[score-minScore] += number; }
	public void addProb(int score, float prob)		{ probDistribution[score-minScore] += prob; }
	
	public float getProbability(int score)
	{
		int index = (score >= minScore) ? score-minScore : 0;
		return probDistribution[index];
	}

	public float getNumberRecs(int score)
	{
		int index = (score >= minScore) ? score-minScore : 0;
		return numDistribution[index];
	}
	
	public float getSpectralProbability(int score)
	{
		float specProb = 0;
		int minIndex = (score >= minScore) ? score-minScore : 0;
		for(int t=minIndex; t<probDistribution.length; t++)
			specProb += probDistribution[t];
		return specProb;
	}

	public float getSpectralProbability(float specProbThreshold)
	{
		float specProb = 0;
		for(int t=probDistribution.length-1; t>=0; t--)
		{
			if(specProb+probDistribution[t] <= specProbThreshold)
				specProb += probDistribution[t];
			else
				break;
		}
		return specProb;
	}
	
	public float getNumEqualOrBetterPeptides(int score)
	{
		float numBetterPeptides = 0;
		int minIndex = (score >= minScore) ? score-minScore : 0;
		for(int t=minIndex; t<numDistribution.length; t++)
			numBetterPeptides += numDistribution[t];
		return numBetterPeptides;
	}
	
	public void addNumDist(ScoreDist otherDist, int scoreDiff)
	{
		addNumDist(otherDist, scoreDiff, 1);
	}

	public void addNumDist(ScoreDist otherDist, int scoreDiff, int coeff)
	{
		if(otherDist == null)
			return;
		for(int t = Math.max(otherDist.minScore, minScore-scoreDiff); t<otherDist.maxScore; t++)
			numDistribution[t+scoreDiff-minScore] += coeff*otherDist.numDistribution[t-otherDist.minScore];
	}
	
	public void addProbDist(ScoreDist otherDist, int scoreDiff, float aaProb)
	{
		if(otherDist == null)
			return;
		for(int t = Math.max(otherDist.minScore, minScore-scoreDiff); t<otherDist.maxScore; t++) {
			float prob = otherDist.probDistribution[t-otherDist.minScore]*aaProb;
			probDistribution[t+scoreDiff-minScore] += prob;	// TODO: underflow
		}
	}
	
	public float getMeanScore()
	{
		double sumScores = 0;
		double sumNum = 0;
		for(int score=this.getMinScore(); score<this.getMaxScore(); score++)
		{
			sumNum += this.getNumberRecs(score);
			sumScores += this.getNumberRecs(score)*score;
		}

		return (float)(sumScores/sumNum);
	}
	
	public ScoreBound getPercentileRange(float percentile)
	{
		return null;
	}
	
	// added by kyowon.  Get a new ScoreDist instance. it has the same value as the original one from newMinScore to max score of the original ScoreDist
	static public ScoreDist getTruncatedScoreDist(ScoreDist original, int newMinScore){
		ScoreDistFactory factory = new ScoreDistFactory(original.isNumSet(), original.isProbSet());
		ScoreDist newDist = factory.getInstance(Math.max(newMinScore, original.getMinScore()), original.getMaxScore());
		
		for(int score = newDist.getMinScore(); score<newDist.getMaxScore(); score++){
			if(newDist.isNumSet()){
				newDist.setNumber(score, original.getNumberRecs(score));
			}
			if(newDist.isProbSet()){
				newDist.setProb(score, original.getProbability(score));
			}
		}
		
		return newDist;
	}
	
}
