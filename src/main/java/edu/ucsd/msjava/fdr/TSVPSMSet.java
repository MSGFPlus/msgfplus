package edu.ucsd.msjava.fdr;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class TSVPSMSet extends PSMSet {

	// required
	File file;
	String delimeter;
	boolean hasHeader;
	int scoreCol;
	boolean isGreaterBetter;
	int specFileCol;
	int specIndexCol;
	int pepCol;
	ArrayList<Pair<Integer,ArrayList<String>>> reqStrList;

	// optional
	int dbCol;
	String decoyPrefix;
	boolean isTarget;

	public TSVPSMSet(
			File file, 
			String delimeter, 
			boolean hasHeader,
			int scoreCol, 
			boolean isGreaterBetter, 
			int specFileCol,
			int specIndexCol, 
			int pepCol,
			ArrayList<Pair<Integer,ArrayList<String>>> reqStrList
	)
	{
		this.file = file;
		this.delimeter = delimeter;
		this.hasHeader = hasHeader;
		this.scoreCol = scoreCol;
		this.isGreaterBetter = isGreaterBetter;
		this.specFileCol = specFileCol;
		this.specIndexCol = specIndexCol;
		this.pepCol = pepCol;
		this.reqStrList = reqStrList;
		dbCol = -1;
	}

	public TSVPSMSet decoy(int dbCol, String decoyPrefix, boolean isTarget)
	{
		this.dbCol = dbCol;
		this.decoyPrefix = decoyPrefix;
		this.isTarget = isTarget;
		return this;
	}

	public String getHeader()	{ return header; }
	public boolean isGreaterBetter()	{ return this.isGreaterBetter; }

	String header;

	public void read()
	{
		psmList = new ArrayList<ScoredString>();		
		peptideScoreTable = new HashMap<String,Float>();

		BufferedReader in = null;
		try {
			in = new BufferedReader(new FileReader(file));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		try {
			if(hasHeader)
			{
				header = in.readLine();
			}
			
			String s;
			HashSet<String> specKeySet = new HashSet<String>();

			while((s=in.readLine()) != null)
			{
				if(s.startsWith("#"))
					continue;
				String[] token = s.split(delimeter);
				if(scoreCol >= token.length || pepCol >= token.length)
					continue;

				String specFile;
				if(specFileCol >= 0)
					specFile = token[specFileCol];
				else
					specFile = "";
				String specIndex = token[specIndexCol];
				String specKey = specFile+":"+specIndex;

				if(specKeySet.contains(specKey))
					continue;
				else
					specKeySet.add(specKey);

				if(dbCol >= 0)
				{
					if(isTarget)
					{
						if(token[dbCol].contains(decoyPrefix))
							continue;
					}
					else
					{
						if(!token[dbCol].contains(decoyPrefix))
							continue;
					}
				}

				if(reqStrList != null)
				{
					boolean isMatched = true;
					for(Pair<Integer,ArrayList<String>> pair : reqStrList)
					{
						boolean containingReqSeq = false;
						for(String reqStr : pair.getSecond())
						{
							if(token[pair.getFirst()].contains(reqStr))
							{
								containingReqSeq = true;
								break;
							}
						}
						if(containingReqSeq == false)
						{
							isMatched = false;
							break;
						}
						else
							isMatched = true;
					}
					if(isMatched == false)
						continue;
				}			

				if(token[scoreCol].length() == 0 || !Character.isDigit(token[scoreCol].charAt(0)))
					continue;
				String pep = getPeptideFromAnnotation(token[pepCol]);
				float score = Float.parseFloat(token[scoreCol]);
				psmList.add(new ScoredString(s, score));

				Float prevScore = peptideScoreTable.get(pep);
				if(prevScore == null || (isGreaterBetter && score > prevScore) || (!isGreaterBetter && score < prevScore))
				{
					peptideScoreTable.put(pep, score);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		if(in != null)
			try {
				in.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
	}

	public void writeResults(TargetDecoyAnalysis tda, PrintStream out, float fdrThreshold, float pepFDRThreshold, float scoreThreshold)
	{
		if(header != null)
			out.println(header + delimeter + "FDR"+delimeter+"PepFDR");
		for(ScoredString ss : getPSMList())
		{
			float psmFDR = tda.getPSMFDR(ss.getScore());
			if(psmFDR > fdrThreshold)
				continue;
			if(isGreaterBetter && ss.getScore() <= scoreThreshold ||
				!isGreaterBetter && ss.getScore() >= scoreThreshold)
				continue;
			String[] token = ss.getStr().split(delimeter);
			Float pepFDR = tda.getPepFDR(token[pepCol]);
			if(pepFDR == null || pepFDR > pepFDRThreshold)
				continue;
			String prevResult = ss.getStr();
			if(!prevResult.endsWith(delimeter))
				prevResult += delimeter;
			out.println(prevResult+psmFDR+delimeter+pepFDR);
		}
		out.flush();
	}
	
	public int getNumIdentifiedPSMs(TargetDecoyAnalysis tda,float fdrThreshold)
	{
		int numID = 0;
		for(ScoredString ss : getPSMList())
		{
			float psmFDR = tda.getPSMFDR(ss.getScore());
			if(psmFDR > fdrThreshold)
				continue;
			numID++;
		}
		return numID;
	}

	public int getNumIdentifiedPeptides(TargetDecoyAnalysis tda, float pepFDRThreshold)
	{
		HashSet<String> pepSet = new HashSet<String>();
		for(ScoredString ss : getPSMList())
		{
			String[] token = ss.getStr().split(delimeter);
			Float pepFDR = tda.getPepFDR(token[pepCol]);
			if(pepFDR == null || pepFDR > pepFDRThreshold)
				continue;
			
			pepSet.add(TSVPSMSet.getPeptideFromAnnotation(token[pepCol]));
		}
		return pepSet.size();
	}
	
	public static String getPeptideFromAnnotation(String annotation)
	{
		String pep;
		// if there are flanking amino acids (e.g. R.ACDEFK.G), remove them
		int firstDotIndex = annotation.indexOf('.');
		int lastDotIndex = annotation.lastIndexOf('.');
		if(firstDotIndex < lastDotIndex)
			pep = annotation.substring(firstDotIndex+1, lastDotIndex);
		else
			pep = annotation;
		pep = pep.toUpperCase();
		return pep;
	}
	
	public static void main(String argv[]) throws Exception
	{
		File file = new File("/home/sangtaekim/Research/ToolDistribution/Test/inspect.out");
		ArrayList<Pair<Integer,ArrayList<String>>> reqStrList = new ArrayList<Pair<Integer,ArrayList<String>>>();
		ArrayList<String> charges = new ArrayList<String>();
		charges.add("1");
		charges.add("3");
		ArrayList<String> peps = new ArrayList<String>();
		peps.add("EE");
		reqStrList.add(new Pair<Integer,ArrayList<String>>(2,peps));
		reqStrList.add(new Pair<Integer,ArrayList<String>>(4,charges));
		TSVPSMSet psmSet = new TSVPSMSet(file, "\t", true, 14, true, 0, 1, 2, reqStrList);
		psmSet.read();
		psmSet.printPeptideScoreTable();
	}	
}
