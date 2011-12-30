package fdr;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.Map.Entry;

public class TargetDecoyPSMSet {
	PSMSet target;
	PSMSet decoy;
	boolean isConcatenated;
	TreeMap<Float,Float> psmLevelFDRMap;	// PSMScore -> FDR
	TreeMap<Float,Float> pepLevelFDRMap;	// Peptide -> PepFDR
	float pit = 1;	// portion of incorrect target PSMs
	
	public TargetDecoyPSMSet(
			File concatenatedFile, 
			String delimeter, 
			boolean hasHeader,
			int scoreCol, 
			boolean isGreaterBetter, 
			int specFileCol,			
			int specIndexCol, 
			int pepCol,
			ArrayList<Pair<Integer,ArrayList<String>>> reqStrList,
			int dbCol, String decoyPrefix)
	{
		target = new PSMSet(concatenatedFile, delimeter, hasHeader, scoreCol, isGreaterBetter, specFileCol, specIndexCol, pepCol, reqStrList).decoy(dbCol, decoyPrefix, true).read();
		decoy = new PSMSet(concatenatedFile, delimeter, hasHeader, scoreCol, isGreaterBetter, specFileCol, specIndexCol, pepCol, reqStrList).decoy(dbCol, decoyPrefix, false).read();
		isConcatenated = true;
		psmLevelFDRMap = getFDRMap(target.getPSMScores(), decoy.getPSMScores(), target.isGreaterBetter, isConcatenated, 1);
		pepLevelFDRMap = getFDRMap(target.getPepScores(), decoy.getPepScores(), target.isGreaterBetter, isConcatenated, 1);
	}
	
	public TargetDecoyPSMSet(
			File targetFile, 
			File decoyFile, 
			String delimeter, 
			boolean hasHeader,
			int scoreCol, 
			boolean isGreaterBetter, 
			int specFileCol,
			int specIndexCol, 
			int pepCol,
			ArrayList<Pair<Integer,ArrayList<String>>> reqStrListPSMSet
			)
	{
		this(targetFile, decoyFile, delimeter, hasHeader, scoreCol, isGreaterBetter, specFileCol, specIndexCol, pepCol, reqStrListPSMSet, 1);
	}
	
	public TargetDecoyPSMSet(
			File targetFile, 
			File decoyFile, 
			String delimeter, 
			boolean hasHeader,
			int scoreCol, 
			boolean isGreaterBetter, 
			int specFileCol,
			int specIndexCol, 
			int pepCol,
			ArrayList<Pair<Integer,ArrayList<String>>> reqStrListPSMSet,
			float pit
			)
	{
		target = new PSMSet(targetFile, delimeter, hasHeader, scoreCol, isGreaterBetter, specFileCol, specIndexCol, pepCol, reqStrListPSMSet).read();
		decoy = new PSMSet(decoyFile, delimeter, hasHeader, scoreCol, isGreaterBetter, specFileCol, specIndexCol, pepCol, reqStrListPSMSet).read();
		isConcatenated = false;
		psmLevelFDRMap = getFDRMap(target.getPSMScores(), decoy.getPSMScores(), target.isGreaterBetter, isConcatenated, pit);		
		pepLevelFDRMap = getFDRMap(target.getPepScores(), decoy.getPepScores(), target.isGreaterBetter, isConcatenated, pit);
	}

	public PSMSet getTargetPSMSet()	{ return target; }
	public PSMSet getDecoyPSMSet() 	{ return decoy; }
	public TreeMap<Float,Float> getPSMLevelFDRMap()	{ return psmLevelFDRMap; }
	public TreeMap<Float,Float> getPepLevelFDRMap()	{ return pepLevelFDRMap; }
	
	public float getPSMFDR(float score)
	{
		float fdr;
		if(target.isGreaterBetter)
			fdr = psmLevelFDRMap.lowerEntry(score).getValue();
		else
			fdr = psmLevelFDRMap.higherEntry(score).getValue();
		return fdr;
	}
	
	public float getPepFDR(float score)
	{
		float fdr;
		if(target.isGreaterBetter)
			fdr = pepLevelFDRMap.lowerEntry(score).getValue();
		else
			fdr = pepLevelFDRMap.higherEntry(score).getValue();
		return fdr;
	}
	
	public Float getPepFDR(String annotation)
	{
		String pep = PSMSet.getPeptideFromAnnotation(annotation);
		
		Float score = target.getPeptideScoreTable().get(pep);
		if(score == null)
		{
			score = decoy.getPeptideScoreTable().get(pep);
			if(score == null)
				return null;
		}
		return getPepFDR(score);
	}
	
	public void writeResults(PrintStream out, float fdrThreshold, float pepFDRThreshold)
	{
		writeResults(out, fdrThreshold, pepFDRThreshold, false);
	}
	
	public void writeResults(PrintStream out, float fdrThreshold, float pepFDRThreshold, boolean showDecoy)
	{
		writeResults(out, fdrThreshold, pepFDRThreshold, (target.isGreaterBetter ? Float.MIN_VALUE : Float.MAX_VALUE), showDecoy);
	}
	
	public void writeResults(PrintStream out, float fdrThreshold, float pepFDRThreshold, float scoreThreshold, boolean showDecoy)
	{
		if(target.getHeader() != null)
			out.println(target.getHeader() + target.delimeter + "FDR"+target.delimeter+"PepFDR");
		for(ScoredString ss : target.getPSMList())
		{
			float psmFDR = getPSMFDR(ss.getScore());
			if(psmFDR > fdrThreshold)
				continue;
			if(target.isGreaterBetter && ss.getScore() <= scoreThreshold ||
				!target.isGreaterBetter && ss.getScore() >= scoreThreshold)
				continue;
			String[] token = ss.getStr().split(target.delimeter);
			Float pepFDR = getPepFDR(token[target.pepCol]);
			if(pepFDR == null || pepFDR > pepFDRThreshold)
				continue;
			String prevResult = ss.getStr();
			if(!prevResult.endsWith(target.delimeter))
				prevResult += target.delimeter;
			out.println(prevResult+psmFDR+target.delimeter+pepFDR);
		}
		if(showDecoy)
		{
			for(ScoredString ss : decoy.getPSMList())
			{
				float psmFDR = getPSMFDR(ss.getScore());
				if(psmFDR > fdrThreshold)
					continue;
				if(decoy.isGreaterBetter && ss.getScore() <= scoreThreshold ||
						!decoy.isGreaterBetter && ss.getScore() >= scoreThreshold)
						continue;
				String[] token = ss.getStr().split(decoy.delimeter);
				Float pepFDR = getPepFDR(token[decoy.pepCol]);
				if(pepFDR == null || pepFDR > pepFDRThreshold)
					continue;
				String prevResult = ss.getStr();
				if(!prevResult.endsWith(decoy.delimeter))
					prevResult += decoy.delimeter;
				out.println(prevResult+psmFDR+decoy.delimeter+pepFDR);
			}
		}
		
		out.flush();
	}

	public int getNumIdentifiedPSMs(float fdrThreshold)
	{
		int numID = 0;
		for(ScoredString ss : target.getPSMList())
		{
			float psmFDR = getPSMFDR(ss.getScore());
			if(psmFDR > fdrThreshold)
				continue;
			numID++;
		}
		return numID;
	}

	// returns threshold where FDR(t>threshold)<=fdrThreshold && FDR(t<=threshold)>fdrThreshold
	public float getThresholdScore(float fdrThreshold, boolean isPeptideLevel)
	{
		TreeMap<Float,Float> map;
		if(!isPeptideLevel)
			map = psmLevelFDRMap;	// PSMScore -> FDR
		else
			map = pepLevelFDRMap;
		
		float threshold;
		if(target.isGreaterBetter)
		{
			threshold = Float.MAX_VALUE;
			for(Entry<Float, Float> entry : map.descendingMap().entrySet())
			{
				if(entry.getValue() > fdrThreshold)
					break;
				else
					threshold = entry.getKey();

			}
		}
		else
		{
			threshold = Float.MIN_VALUE;
			
			for(Entry<Float, Float> entry : map.entrySet())
			{
//				System.out.println(entry.getKey()+"\t"+entry.getValue());
				if(entry.getValue() > fdrThreshold)
					break;
				else
					threshold = entry.getKey();
			}
		}
		return threshold;
	}
	
	public int getNumIdentifiedPeptides(float pepFDRThreshold)
	{
		HashSet<String> pepSet = new HashSet<String>();
		for(ScoredString ss : target.getPSMList())
		{
			String[] token = ss.getStr().split(target.delimeter);
			Float pepFDR = getPepFDR(token[target.pepCol]);
			if(pepFDR == null || pepFDR > pepFDRThreshold)
				continue;
			
			pepSet.add(PSMSet.getPeptideFromAnnotation(token[target.pepCol]));
		}
		return pepSet.size();
	}
	
	public static TreeMap<Float,Float> getFDRMap(ArrayList<Float> target, ArrayList<Float> decoy, 
			boolean isGreaterBetter, boolean isConcatenated, float pit)
	{
		TreeMap<Float,Float> fdrMap = new TreeMap<Float,Float>();
		if(!isGreaterBetter)
		{
			Collections.sort(target);
			Collections.sort(decoy);
		}
		else
		{
			Collections.sort(target, Collections.reverseOrder());
			Collections.sort(decoy, Collections.reverseOrder());
		}
		
		int targetIndex = 0;
		float prevDecoyScore = Float.MIN_VALUE;
		
		for(int decoyIndex=0; decoyIndex<decoy.size(); decoyIndex++)
		{
			float decoyScore = decoy.get(decoyIndex);
			if(decoyScore == prevDecoyScore)
				continue;
			else
				prevDecoyScore = decoyScore;
			if(isGreaterBetter)
			{
				while(targetIndex < target.size() && target.get(targetIndex) > decoyScore)
					targetIndex++;
			}
			else
			{
				while(targetIndex < target.size() && target.get(targetIndex) < decoyScore)
					targetIndex++;
			}
			
			if(targetIndex > 0)
			{
				float fdr;
				if(targetIndex <= decoyIndex)
					fdr = 1;
				else
				{
					if(!isConcatenated)
					{
						fdr = Math.round(decoyIndex*pit)/(float)targetIndex;	// simple formulation by Kall et et., JPR 2008
//						fdr = (2*decoyIndex)/(float)(targetIndex + decoyIndex);	// Elias and Gygi, Nat. Methods 2007
					}
					else
//						fdr = (2*decoyIndex)/(float)(targetIndex + decoyIndex);	// Elias and Gygi, Nat. Methods 2007
						fdr = decoyIndex/(float)targetIndex;
				}
				if(fdr > 1)
					fdr = 1f;

				fdrMap.put(decoyScore, fdr);
				if(fdr >= 1)
					break;
//				System.out.println("1: " + decoyScore+":"+fdr);
			}
		}
		
		TreeMap<Float,Float> finalFDRMap = new TreeMap<Float,Float>();
		
		// Convert FDRs into q-values
		Iterator<Entry<Float, Float>> itr;
		if(isGreaterBetter)
			itr = fdrMap.entrySet().iterator();
		else
			itr = fdrMap.descendingMap().entrySet().iterator();
		float minFDR = 1;
		while(itr.hasNext())
		{
			Entry<Float,Float> entry = itr.next();
			float fdr = entry.getValue();
			if(fdr > minFDR)
				fdr = minFDR;
			minFDR = fdr;
			finalFDRMap.put(entry.getKey(), fdr);
		}
		if(isGreaterBetter)
			finalFDRMap.put(Float.NEGATIVE_INFINITY, 1f);
		else
			finalFDRMap.put(Float.POSITIVE_INFINITY, 1f);
		return finalFDRMap;
	}	
	
}
