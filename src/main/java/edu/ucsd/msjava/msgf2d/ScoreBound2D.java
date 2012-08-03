package edu.ucsd.msjava.msgf2d;

import edu.ucsd.msjava.msgf.ScoreBound;

public class ScoreBound2D {
	protected ScoreBound scoreBound1;
	protected ScoreBound scoreBound2;
	
	public ScoreBound2D(int minScore1, int maxScore1, int minScore2, int maxScore2)
	{
		this.scoreBound1 = new ScoreBound(minScore1, maxScore1);
		this.scoreBound2 = new ScoreBound(minScore2, maxScore2);
	}
	public ScoreBound2D(ScoreBound2D scoreBound) {
		this.scoreBound1 = scoreBound.scoreBound1;
		this.scoreBound2 = scoreBound.scoreBound2;
	}
	public ScoreBound2D(ScoreBound scoreBound1, ScoreBound scoreBound2) {
		this.scoreBound1 = scoreBound1;
		this.scoreBound2 = scoreBound2;
	}
	public int getMinScore1() {
		return scoreBound1.getMinScore();
	}
	public int getMinScore2() {
		return scoreBound2.getMinScore();
	}
	public void setMinScore(int minScore1, int minScore2) {
		scoreBound1.setMinScore(minScore1);
		scoreBound2.setMinScore(minScore2);
	}
	public void setMaxScore(int maxScore1, int maxScore2) {
		scoreBound1.setMaxScore(maxScore1);
		scoreBound2.setMaxScore(maxScore2);
	}
	public int getMaxScore1() {
		return scoreBound1.getMaxScore();
	}
	public int getMaxScore2() {
		return scoreBound2.getMaxScore();
	}
	public int getRange1() {
		return getMaxScore1() - getMinScore1();
	}
	public int getRange2() {
		return getMaxScore2() - getMinScore2();
	}

}
