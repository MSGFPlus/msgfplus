package msgap;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;

import msgf.ScoreBound;
import msutil.Matter;

public class GapBacktrackPointer<T extends Matter>  extends ScoreBound{
	
	ArrayList<T> hubSet = null;
	ArrayList<Integer>[] connectedHubIndicesTable = null;
	private BitSet[] backtrackPointer = null; 
	private int[] scoreDiffs = null;
	
	@SuppressWarnings("unchecked")
	public GapBacktrackPointer(ArrayList<T> hubSet, int minScore, int maxScore){
		super(minScore, maxScore);
		this.hubSet = hubSet;
		scoreDiffs = new int[hubSet.size()];
		
		connectedHubIndicesTable = new ArrayList[maxScore - minScore];
		for(int i=0;i<maxScore - minScore; i++)
			connectedHubIndicesTable[i] = new ArrayList<Integer>();
		for(int i=0; i<scoreDiffs.length; i++) scoreDiffs[i] = Integer.MIN_VALUE;
		backtrackPointer = new BitSet[maxScore - minScore];
		for(int i=0; i<backtrackPointer.length; i++) backtrackPointer[i] = new BitSet(hubSet.size());
	}
	
	public int getScoreDiff(int hubIndex) {
		return scoreDiffs[hubIndex]; 
	}

	public ArrayList<Integer> getConnectedHubIndices(int score) {
		return connectedHubIndicesTable[score-minScore];
	}
	
	public int getScoreDiff(T hub) {
		return getScoreDiff(hubSet.indexOf(hub)); 
	}
	
	public boolean isSet(int score, int scoreDiff, int hubIndex)
	{
		if(scoreDiff == Integer.MIN_VALUE) return false;
		if(score < minScore || score >= maxScore) return false;
		if(scoreDiffs[hubIndex] != scoreDiff) return false;
		return backtrackPointer[score-minScore].get(hubIndex);
	}
	
	/*public boolean isSet(int score, int scoreDiff, T hub)
	{
		return this.isSet(score, scoreDiff, hubSet.indexOf(hub));
	}*/
	
	public void clearBacktrack(int score, int hubIndex) { 
		backtrackPointer[score-minScore].clear(hubIndex); 
	}
	public void clearBacktrack(int score, T hub) { 
		backtrackPointer[score-minScore].clear(hubSet.indexOf(hub));
	}
	
	public void setBacktrack(int score, int scoreDiff, int hubIndex)	{ 
		backtrackPointer[score-minScore].set(hubIndex);
		scoreDiffs[hubIndex] = scoreDiff;
		if(!connectedHubIndicesTable[score-minScore].contains(hubIndex))
			connectedHubIndicesTable[score-minScore].add(hubIndex);
	}
	public void setBacktrack(int score, int scoreDiff, T hub)	{
		setBacktrack(score, scoreDiff, hubSet.indexOf(hub));
	}
	
	public BitSet getBacktrackPointers(int score)	{ return backtrackPointer[score-minScore]; }
	
	public void addBacktrackPointers(GapBacktrackPointer<T> prevPointer, int hubIndex, int scorediff)
	{
		for(int t = Math.max(prevPointer.minScore, minScore-scorediff); t<prevPointer.maxScore; t++)
		{
			if(!prevPointer.getBacktrackPointers(t).isEmpty()){
				this.setBacktrack(t+scorediff, scorediff, hubIndex);
			}
		}
	}
}
