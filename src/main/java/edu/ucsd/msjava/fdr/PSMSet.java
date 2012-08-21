package edu.ucsd.msjava.fdr;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

public abstract class PSMSet {
	protected ArrayList<ScoredString> psmList;	// resultLine, psm
	protected HashMap<String,Float> peptideScoreTable;	// peptide -> best score

	public ArrayList<ScoredString> getPSMList()	{ return psmList; }
	public HashMap<String,Float> getPeptideScoreTable() { return peptideScoreTable; }
	public abstract boolean isGreaterBetter();
	
	public void printPSMSet()
	{
		if(psmList != null)
		{
			for(ScoredString s : psmList)
			{
				System.out.println(s.getStr());
			}
		}
	}

	public void printPeptideScoreTable()
	{
		if(peptideScoreTable != null)
		{
			Iterator<Entry<String, Float>> itr = peptideScoreTable.entrySet().iterator();
			while(itr.hasNext())
			{
				Entry<String,Float> entry = itr.next();
				System.out.println(entry.getKey()+"\t"+entry.getValue());
			}
		}
	}

	public ArrayList<Float> getPSMScores()
	{
		if(psmList == null)
			return null;
		ArrayList<Float> psmScores = new ArrayList<Float>();
		for(ScoredString ss : psmList)
			psmScores.add(ss.getScore());
		return psmScores;
	}
	
	public ArrayList<Float> getPepScores()
	{
		if(peptideScoreTable == null)
			return null;
		ArrayList<Float> pepScores = new ArrayList<Float>();
		Iterator<Entry<String, Float>> itr = peptideScoreTable.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<String,Float> entry = itr.next();
			pepScores.add(entry.getValue());
		}
		return pepScores;
	}
	
	public void writeResults(TargetDecoyAnalysis tda, PrintStream out, float fdrThreshold, float pepFDRThreshold)
	{
		writeResults(tda, out, fdrThreshold, pepFDRThreshold);
	}
	
	public abstract void writeResults(TargetDecoyAnalysis tda, PrintStream out, float fdrThreshold, float pepFDRThreshold, float scoreThreshold);
}
