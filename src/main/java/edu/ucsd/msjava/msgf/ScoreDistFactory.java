package edu.ucsd.msjava.msgf;

public class ScoreDistFactory {
	boolean calcNumber, calcProb;
	
	public ScoreDistFactory(boolean calcNumber, boolean calcProb) 
	{
		this.calcNumber = calcNumber;
		this.calcProb = calcProb;
	}
	
	public ScoreDist getInstance(int minScore, int maxScore)
	{
		return new ScoreDist(minScore, maxScore, calcNumber, calcProb);
	}
}
