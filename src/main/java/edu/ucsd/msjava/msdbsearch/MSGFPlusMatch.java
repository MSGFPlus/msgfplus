package edu.ucsd.msjava.msdbsearch;

import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;

public class MSGFPlusMatch implements Comparable<MSGFPlusMatch> {

	private final int specIndex;
	private final List<DatabaseMatch> matchList;
	private final double specEValue;
	
	public MSGFPlusMatch(int specIndex, PriorityQueue<DatabaseMatch> matchQueue)
	{
		this.specIndex = specIndex;
		this.matchList = new ArrayList<DatabaseMatch>(matchQueue);
		specEValue = matchList.get(matchList.size()-1).getSpecEValue();
	}

	public int getSpecIndex() {
		return specIndex;
	}

	public List<DatabaseMatch> getMatchList() {
		return matchList;
	}

	public double getSpecEValue() {
		return specEValue;
	}
	
	@Override
	public int compareTo(MSGFPlusMatch o) {
		if(specEValue < o.specEValue)
			return -1;
		else if(specEValue == o.specEValue)
			return 0;
		else
			return 1;
	}
	
}
