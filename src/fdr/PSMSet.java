package fdr;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map.Entry;

public class PSMSet {

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

	public PSMSet(File file, String delimeter, boolean hasHeader,
			int scoreCol, boolean isGreaterBetter, 
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

	public PSMSet decoy(int dbCol, String decoyPrefix, boolean isTarget)
	{
		this.dbCol = dbCol;
		this.decoyPrefix = decoyPrefix;
		this.isTarget = isTarget;
		return this;
	}

	public String getHeader()	{ return header; }
	public ArrayList<ScoredString> getPSMList()	{ return psmList; }
	public HashMap<String,Float> getPeptideScoreTable() { return peptideScoreTable; }
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

	String header;
	ArrayList<ScoredString> psmList;	// resultLine, psm
	HashMap<String,Float> peptideScoreTable;	// peptide -> best score

	public PSMSet read()
	{
		psmList = new ArrayList<ScoredString>();		
		peptideScoreTable = new HashMap<String,Float>();

		int lineNum = 0;
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
				lineNum++;
			}
			
			String s;
			HashSet<String> specKeySet = new HashSet<String>();

			while((s=in.readLine()) != null)
			{
				lineNum++;
				
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
				
				///////////////////////////////////////////
//				if(pep.contains("KKKK") || pep.contains("DDDDDDD"))
//					continue;
//				int numMissedCleavages = 0;
//				for(int i=0; i<pep.length(); i++)
//				{
//					if(pep.charAt(i) == 'K')
//						numMissedCleavages++;
//				}
////				if(numMissedCleavages > Math.round(pep.length()*0.3f))
////					continue;
//				if(numMissedCleavages > Math.round(pep.length()*0.3f))
//					continue;
//				// debug
//				float pmError = Float.parseFloat(token[2]);
//				int charge = Integer.parseInt(token[3]);
//				if(pmError > 5)
//					continue;
				//////////////////////////////////////
				
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
		return this;
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
		PSMSet psmSet = new PSMSet(file, "\t", true, 14, true, 0, 1, 2, reqStrList);
		psmSet.read();
		psmSet.printPeptideScoreTable();
	}	
}
