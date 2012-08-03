package edu.ucsd.msjava.msgf;

public class ScoreBound {
	protected int minScore;	// inclusive
	protected int maxScore;	// exclusive
	public ScoreBound(int minScore, int maxScore) {
		this.minScore = minScore;
		this.maxScore = maxScore;
	}
	public int getMinScore() {
		return minScore;
	}
	public void setMinScore(int minScore) {
		this.minScore = minScore;
	}
	public int getMaxScore() {
		return maxScore;
	}
	public int getRange() {
		return maxScore-minScore;
	}
	public void setMaxScore(int maxScore) {
		this.maxScore = maxScore;
	}

}
