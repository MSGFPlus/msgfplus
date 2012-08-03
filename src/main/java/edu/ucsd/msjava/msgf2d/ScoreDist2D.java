package edu.ucsd.msjava.msgf2d;

public class ScoreDist2D extends ScoreBound2D {
	private float[][] numDistribution;
	private float[][] probDistribution;
	
	public ScoreDist2D(int minScore1, int maxScore1, int minScore2, int maxScore2)
	{
		super(minScore1, maxScore1, minScore2, maxScore2);
		numDistribution = new float[scoreBound1.getRange()][scoreBound2.getRange()];
		probDistribution = new float[scoreBound1.getRange()][scoreBound2.getRange()];
	}
	
	public ScoreDist2D(ScoreBound2D scoreBound)
	{
		super(scoreBound.scoreBound1, scoreBound.scoreBound2);
		numDistribution = new float[scoreBound1.getRange()][scoreBound2.getRange()];
		probDistribution = new float[scoreBound1.getRange()][scoreBound2.getRange()];
	}
	
	public void setNumber(int score1, int score2, float number)	{ numDistribution[score1-scoreBound1.getMinScore()][score2-scoreBound2.getMinScore()] = number; }
	public void setProb(int score1, int score2, float prob)	{ probDistribution[score1-scoreBound1.getMinScore()][score2-scoreBound2.getMinScore()] = prob; }
	
	public void addNumber(int score1, int score2, float number)	{ numDistribution[score1-scoreBound1.getMinScore()][score2-scoreBound2.getMinScore()] += number; }
	public void addProb(int score1, int score2, float prob)	{ probDistribution[score1-scoreBound1.getMinScore()][score2-scoreBound2.getMinScore()] += prob; }

	public float getProbability(int score1, int score2)
	{
		int index1 = (score1 >= scoreBound1.getMinScore()) ? score1-scoreBound1.getMinScore() : 0;
		int index2 = (score2 >= scoreBound2.getMinScore()) ? score2-scoreBound2.getMinScore() : 0;
		return probDistribution[index1][index2];
	}
	
	public float getNumRecs(int score1, int score2)
	{
		int index1 = (score1 >= scoreBound1.getMinScore()) ? score1-scoreBound1.getMinScore() : 0;
		int index2 = (score2 >= scoreBound2.getMinScore()) ? score2-scoreBound2.getMinScore() : 0;
		return numDistribution[index1][index2];
	}
		
	/**
	 * @deprecated
	 */
	public float getSpectralProbability1(int score1)
	{
		int index1 = (score1 >= scoreBound1.getMinScore()) ? score1-scoreBound1.getMinScore() : 0;
		float specProb = 0;
		for(int t1=index1; t1<scoreBound1.getRange(); t1++)
		{
			for(int t2=0; t2<scoreBound2.getRange(); t2++)
			{
				specProb += probDistribution[t1][t2];
			}
		}
			
		return specProb;
	}

	/**
	 * @deprecated
	 */
	public float getSpectralProbability2(int score2)
	{
		int index2 = (score2 >= scoreBound2.getMinScore()) ? score2-scoreBound2.getMinScore() : 0;
		float specProb = 0;
		for(int t1=0; t1<scoreBound1.getRange(); t1++)
		{
			for(int t2=index2; t2<scoreBound2.getRange(); t2++)
			{
				specProb += probDistribution[t1][t2];
			}
		}
		return specProb;
	}
	
	/**
	 * @deprecated
	 */
	public float getSpectralProbabilitySumScores(int score1, int score2)
	{
		float specProb = 0;
		for(int s1=scoreBound1.getMinScore(); s1<scoreBound1.getMaxScore(); s1++)
		{
			for(int s2=score1+score2-s1; s2<scoreBound2.getMaxScore(); s2++)
			{
				if(s1<scoreBound1.getMinScore() || s2<scoreBound2.getMinScore())
					continue;
				if(s1 >= scoreBound1.getMinScore() && s1 < scoreBound1.getMaxScore() &&
						s2 >= scoreBound2.getMinScore() && s2 < scoreBound2.getMaxScore())
					specProb += probDistribution[s1-scoreBound1.getMinScore()][s2-scoreBound2.getMinScore()];
			}
		}
		return specProb;
	}
	
	/**
	 * @deprecated
	 */
	public float getSpectralProbability(int score1, int score2)
	{
		int index1 = (score1 >= scoreBound1.getMinScore()) ? score1-scoreBound1.getMinScore() : 0;
		int index2 = (score2 >= scoreBound2.getMinScore()) ? score2-scoreBound2.getMinScore() : 0;
		float specProb = 0;
		for(int t1=index1; t1<scoreBound1.getRange(); t1++)
		{
			for(int t2=index2; t2<scoreBound2.getRange(); t2++)
			{
				specProb += probDistribution[t1][t2];
			}
		}
			
		return specProb;
	}

	/**
	 * @deprecated
	 */
	public float getNumEqualOrBetterPeptides(int score1, int score2)
	{
		int index1 = (score1 >= scoreBound1.getMinScore()) ? score1-scoreBound1.getMinScore() : 0;
		int index2 = (score2 >= scoreBound2.getMinScore()) ? score2-scoreBound2.getMinScore() : 0;
		float numRecs = 0;
		for(int t1=index1; t1<scoreBound1.getRange(); t1++)
			for(int t2=index2; t2<scoreBound2.getRange(); t2++)
				numRecs += numDistribution[t1][t2];
			
		return numRecs;
	}
	
	/**
	 * @deprecated
	 */
	public float getNumEqualOrBetterPeptidesSumScores(int score1, int score2)
	{
		float numRecs = 0;
		for(int s1=scoreBound1.getMinScore(); s1<scoreBound1.getMaxScore(); s1++)
		{
			for(int s2=score1+score2-s1; s2<scoreBound2.getMaxScore(); s2++)
			{
				if(s1 >= scoreBound1.getMinScore() && s1 < scoreBound1.getMaxScore() &&
						s2 >= scoreBound2.getMinScore() && s2 < scoreBound2.getMaxScore())
					numRecs += numDistribution[s1-scoreBound1.getMinScore()][s2-scoreBound2.getMinScore()];
			}
		}
		return numRecs;
	}	
	
	/**
	 * @deprecated
	 */
	public float getNumEqualOrBetterPeptides1(int score1)
	{
		int index1 = (score1 >= scoreBound1.getMinScore()) ? score1-scoreBound1.getMinScore() : 0;
		float numRecs = 0;
		for(int t1=index1; t1<scoreBound1.getRange(); t1++)
		{
			for(int t2=0; t2<scoreBound2.getRange(); t2++)
			{
				numRecs += numDistribution[t1][t2];
//				System.out.println("***"+t1+"\t"+t2+"\t"+numDistribution[t1][t2]);
			}
		}
		return numRecs;
	}

	/**
	 * @deprecated
	 */
	public float getNumEqualOrBetterPeptides2(int score2)
	{
		int index2 = (score2 >= scoreBound2.getMinScore()) ? score2-scoreBound2.getMinScore() : 0;
		float numRecs = 0;
		for(int t1=0; t1<scoreBound1.getRange(); t1++)
		{
			for(int t2=index2; t2<scoreBound2.getRange(); t2++)
			{
				numRecs += numDistribution[t1][t2];
//				System.out.println("***"+t1+"\t"+t2+"\t"+numDistribution[t1][t2]);
			}
		}
		return numRecs;
	}
	
	public void addNumDist(ScoreDist2D otherDist, int scoreDiff1, int scoreDiff2, int coeff)
	{
		if(otherDist == null)
			return;
		int minScore1 = scoreBound1.getMinScore();
		int minScore2 = scoreBound2.getMinScore();
		
		for(int t1 = Math.max(otherDist.scoreBound1.getMinScore(), scoreBound1.getMinScore()-scoreDiff1); t1<otherDist.scoreBound1.getMaxScore(); t1++)
		{
			for(int t2 = Math.max(otherDist.scoreBound2.getMinScore(), scoreBound2.getMinScore()-scoreDiff2); t2<otherDist.scoreBound2.getMaxScore(); t2++)
			{
				numDistribution[t1+scoreDiff1-minScore1][t2+scoreDiff2-minScore2] += coeff*otherDist.numDistribution[t1-otherDist.scoreBound1.getMinScore()][t2-otherDist.scoreBound2.getMinScore()];
			}
		}
	}
	public void addProbDist(ScoreDist2D otherDist, int scoreDiff1, int scoreDiff2, float aaProb)
	{
		if(otherDist == null)
			return;
		int minScore1 = scoreBound1.getMinScore();
		int minScore2 = scoreBound2.getMinScore();
		
		for(int t1 = Math.max(otherDist.scoreBound1.getMinScore(), scoreBound1.getMinScore()-scoreDiff1); t1<otherDist.scoreBound1.getMaxScore(); t1++)
		{
			for(int t2 = Math.max(otherDist.scoreBound2.getMinScore(), scoreBound2.getMinScore()-scoreDiff2); t2<otherDist.scoreBound2.getMaxScore(); t2++)
			{
				probDistribution[t1+scoreDiff1-minScore1][t2+scoreDiff2-minScore2] += aaProb*otherDist.probDistribution[t1-otherDist.scoreBound1.getMinScore()][t2-otherDist.scoreBound2.getMinScore()];
			}
		}
	}

}
