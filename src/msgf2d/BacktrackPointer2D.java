package msgf2d;

public class BacktrackPointer2D extends ScoreBound2D {
	private int[][] backtrackPointer;
	private int curScore1;
	private int curScore2;
	
	// minScore: inclusive, maxScore: exclusive
	public BacktrackPointer2D(ScoreBound2D scoreBound, int curScore1, int curScore2)
	{
		super(scoreBound);
		this.curScore1 = curScore1;
		this.curScore2 = curScore2;
		backtrackPointer = new int[scoreBound1.getRange()][scoreBound2.getRange()];
	}
	
	public int getCurScore1()	{ return curScore1; }
	public int getCurScore2()	{ return curScore2; }
	
	public void setBacktrack(int score1, int score2, int aaIndex)	{ backtrackPointer[score1-scoreBound1.getMinScore()][score2-scoreBound2.getMinScore()] |= (1 << aaIndex);}
	public int getBacktrackPointers(int score1, int score2)	{ return backtrackPointer[score1-scoreBound1.getMinScore()][score2-scoreBound2.getMinScore()]; }
	
	public boolean isSet(int score1, int score2, int aaIndex)
	{
		int mask = (1 << aaIndex);
		if((backtrackPointer[score1-scoreBound1.getMinScore()][score2-scoreBound2.getMinScore()] & mask) == 0)
			return false;
		else
			return true;
	}
	
	public void addBacktrackPointers(BacktrackPointer2D prevPointer, int aaIndex)
	{
		for(int t1 = Math.max(prevPointer.scoreBound1.getMinScore(), scoreBound1.getMinScore()-curScore1); t1<prevPointer.scoreBound1.getMaxScore(); t1++)
		{
			for(int t2 = Math.max(prevPointer.scoreBound2.getMinScore(), scoreBound2.getMinScore()-curScore2); t2<prevPointer.scoreBound2.getMaxScore(); t2++)
			{
				if(prevPointer.getBacktrackPointers(t1, t2) != 0)
					this.setBacktrack(t1+curScore1, t2+curScore2, aaIndex);
			}
		}
	}
}

