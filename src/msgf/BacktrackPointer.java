package msgf;

import java.util.ArrayList;

public class BacktrackPointer extends ScoreBound {
	private int[] backtrackPointer;
	int nodeScore;
	// minScore: inclusive, maxScore: exclusive
	public BacktrackPointer(int minScore, int maxScore, int curScore)
	{
		super(minScore, maxScore);
		this.nodeScore = curScore;
		backtrackPointer = new int[maxScore - minScore];
	}
	
	public int getNodeScore()	{ return nodeScore; }
	public void setBacktrack(int score, int aaIndex)	{ backtrackPointer[score-minScore] |= (1 << aaIndex);}
	
	public int getBacktrackPointers(int score)	{ return backtrackPointer[score-minScore]; }
	
	public boolean isSet(int score, int aaIndex)
	{
		int mask = (1 << aaIndex);
		if((backtrackPointer[score-minScore] & mask) == 0)
			return false;
		else
			return true;
	}
	
	public void addBacktrackPointers(BacktrackPointer prevPointer, int aaIndex, int edgeScore)
	{
		int combinedScore = nodeScore+edgeScore;
		for(int t = Math.max(prevPointer.minScore, minScore-combinedScore); t<prevPointer.maxScore; t++)
		{
			if(prevPointer.getBacktrackPointers(t) != 0)
				this.setBacktrack(t+combinedScore, aaIndex);
		}
	}
	
	public ArrayList<Integer> getBacktrackAAIndexList(int score)
	{
		assert(score >= minScore && score < maxScore);
		int pointer = backtrackPointer[score-minScore];
		int mask = 1;
		ArrayList<Integer> prevIndexList = new ArrayList<Integer>();

		for(int i=0; pointer != 0; i++)
		{
			//if((pointer & (mask << i) ) != 0)
			if((pointer & mask) != 0)
				prevIndexList.add(i);
				
			pointer = pointer >>> 1;
		}
		return prevIndexList;
	}
}
