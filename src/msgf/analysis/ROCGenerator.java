package msgf.analysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.Map.Entry;

import msutil.Pair;
import msutil.ScoredString;
import parser.BufferedLineReader;

public class ROCGenerator {
	public static final float FDR_REPORT_THRESHOLD = 0.1f;
	public static void main(String argv[])
	{
		// required
		File targetFile = null;
		int scoreCol = -1;
		String identifier = null;
		
		// optional
		File outputFile = null;
		boolean isGreaterBetter = false;
		boolean hasHeader = true;
		File decoyFile = null;
		String delimeter = "\t";
		int pepCol = -1;
		int scanNumCol = -1;
		boolean isConcatenated = false;
		
		int dbCol = -1;
		String decoyPrefix = null;
		
		ArrayList<Pair<Integer,String>> reqStrList = null;
		
		int i=0;
		while(i<argv.length)
		{
     		if(argv[i].equalsIgnoreCase("-i"))
     		{
     			identifier = argv[i+1];
     			i+=2;
     		}
     		// 	-f resuleFileName dbCol decoyPrefix or -f targetFileName decoyFileName 
     		else if(argv[i].equalsIgnoreCase("-f"))
			{
				if(i+2 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				targetFile = new File(argv[i+1]);
				if(!targetFile.exists())
					printUsageAndExit(argv[i+1] + " doesn't exist.");
				else if(!targetFile.isFile())
					printUsageAndExit(argv[i+1] + " is not a file.");
				if(i+3 < argv.length && !argv[i+3].startsWith("-"))	// concatenated; -f resultFileName dbCol decoyPrefix
				{
					dbCol = Integer.parseInt(argv[i+2]);
					decoyPrefix = argv[i+3];
					isConcatenated = true;
					i+=4;
				}
				else	// separate; -f targetFileName decoyFileName
				{
					decoyFile = new File(argv[i+2]);
					if(!decoyFile.exists())
						printUsageAndExit(argv[i+2] + " doesn't exist.");
					else if(!decoyFile.isFile())
						printUsageAndExit(argv[i+2] + " is not a file.");
					isConcatenated = false;
					i+=3;
				}
			}
			else if(argv[i].equalsIgnoreCase("-s"))
			{
				if(i+2 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				
				try {
					scoreCol = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal scoreCol: " + argv[i+1]);
				}
				if(argv[i+2].equalsIgnoreCase("1"))
					isGreaterBetter = true;
				else
					isGreaterBetter = false;
				i+=3;
			}
			else if(argv[i].equalsIgnoreCase("-o"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				outputFile = new File(argv[i+1]);
			}
			else if(argv[i].equalsIgnoreCase("-h"))
			{
				if(argv[i+1].equalsIgnoreCase("0"))
					hasHeader = false;
				i+=2;
			}
			else if(argv[i].equalsIgnoreCase("-delim"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				delimeter = argv[i+1];
				i+=2;
			}
			else if(argv[i].equalsIgnoreCase("-p"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				try {
					pepCol = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal pepCol: " + argv[i+1]);
				}
				i+=2;
			}
			else if(argv[i].equalsIgnoreCase("-n"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				try {
					scanNumCol = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal pepCol: " + argv[i+1]);
				}
				i+=2;
			}
			else if(argv[i].equalsIgnoreCase("-m"))
			{
				if(reqStrList == null)
					reqStrList = new ArrayList<Pair<Integer,String>>();
				int matchCol = -1;
				if(i+2 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				try {
					matchCol = Integer.parseInt(argv[i+1]);
				} catch (NumberFormatException e)
				{
					printUsageAndExit("Illigal matchCol: " + argv[i+1]);
				}
				String[] token = argv[i+2].split(",");
				for(String requiredStr : token)
					reqStrList.add(new Pair<Integer,String>(matchCol,requiredStr));
				i+=3;
			}
			else
			{
				printUsageAndExit("Illegal parameter");
			}
		}
		
		if(targetFile == null)
			printUsageAndExit("Target is missing!");
		if(scoreCol < 0)
			printUsageAndExit("scoreCol is missing or illegal!");
		
		printROCCurve(identifier, targetFile, decoyFile, 
				scoreCol, isGreaterBetter, 
				delimeter, scanNumCol, pepCol, reqStrList, 
				isConcatenated, hasHeader, dbCol, decoyPrefix, outputFile);
	}
	
	public static void printUsageAndExit(String message)
	{
		System.err.println(message);
		System.out.print("usage: java ROCGenerator \n" +
				"\t -f resuleFileName dbCol decoyPrefix or -f targetFileName decoyFileName\n" +
				"\t -s scoreCol 0/1 (0: smaller better, 1: greater better)\n" +
				"\t [-o outputFile]\n" +
				"\t [-delim delimeter] (default: \\t)\n" +
				"\t [-p pepCol] (if specified, the peptide level FDRs will be calculated)\n" +
				"\t [-n scanNumCol] (if specified, only best score per spectrum will be considered)\n" +
				"\t [-m colNum keyword (the column 'colNum' must contain 'keyword'. If 'keyword' is delimetered by '|' (e.g. A,B,C), then at least one must be matched.)]\n" +
				"\t [-h 0/1] (0: no header, 1: header (default))\n" +
				"\t [-i identifier (to generate a Matlab code for ROC curves]\n"
				);
		System.exit(-1);
	}

	public static void printROCCurve(String identifier, File targetFile, File decoyFile, int scoreCol, boolean isGreaterBetter, String delimeter, 
			int scanNumCol, int pepCol, ArrayList<Pair<Integer,String>> reqStrList, boolean isConcatenated, boolean hasHeader, int dbCol, String decoyPrefix,
			File outputFile)
	{
		ArrayList<Float> target = null;
		ArrayList<Float> decoy = null;
		
		if(dbCol >= 0)	// both target and decoy are in the same file
		{
			target = getScoreList(getScoredStringList(targetFile, scoreCol, isGreaterBetter, delimeter, scanNumCol, pepCol, reqStrList, hasHeader, dbCol, decoyPrefix, true));
			decoy = getScoreList(getScoredStringList(targetFile, scoreCol, isGreaterBetter, delimeter, scanNumCol, pepCol, reqStrList, hasHeader, dbCol, decoyPrefix, false));
		}
		else
		{
			target = getScoreList(getScoredStringList(targetFile, scoreCol, isGreaterBetter, delimeter, scanNumCol, pepCol, reqStrList, hasHeader, dbCol, decoyPrefix, true));
			decoy = getScoreList(getScoredStringList(decoyFile, scoreCol, isGreaterBetter, delimeter, scanNumCol, pepCol, reqStrList, hasHeader, dbCol, decoyPrefix, true));
		}
		printROCCurve(identifier, target, decoy, isGreaterBetter, isConcatenated);
	}

	private static ArrayList<Float> getScoreList(ArrayList<ScoredString> list)
	{
		ArrayList<Float> scoreList = new ArrayList<Float>();
		for(ScoredString s : list)
			scoreList.add(s.getScore());
		return scoreList;
	}
	
	private static ArrayList<Float> getScoreList(Hashtable<String, Float> table)
	{
		ArrayList<Float> list = new ArrayList<Float>();
		Iterator<Entry<String, Float>> itr = table.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<String, Float> entry = itr.next();
			list.add(entry.getValue());
		}
		return list;
	}
	
	private static ArrayList<ScoredString> getScoredStringList(File file, int scoreCol, boolean isGreaterBetter, String delimeter, int scanNumCol, 
			int pepCol, ArrayList<Pair<Integer,String>> reqStrList, boolean hasHeader, int dbCol, String decoyPrefix, boolean isTarget)
	{
		ArrayList<ScoredString> list = new ArrayList<ScoredString>();
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(file.getPath());
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		if(hasHeader)
			in.readLine();
		
		Hashtable<String, Float> prevScoreTable = new Hashtable<String, Float>();
		ArrayList<Float> scoreList = new ArrayList<Float>();
		String s;
		String prevScanNum = "asdfasfdasdf";	// randomString
		
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#"))
				continue;
			String[] token = s.split(delimeter);
			if(scoreCol >= token.length || pepCol >= token.length || dbCol >= token.length)
				continue;
			
			if(scanNumCol >= 0)
			{
				String scanNum = token[scanNumCol];
				
				if(scanNum.equalsIgnoreCase(prevScanNum))
					continue;
				else	// new scan
				{
					prevScanNum = scanNum;
				}
			}
			if(dbCol >= 0)
			{
				if(isTarget)
				{
					if(token[dbCol].startsWith(decoyPrefix))
						continue;
				}
				else	// decoy
				{
					if(!token[dbCol].startsWith(decoyPrefix))
						continue;
				}
			}
			if(reqStrList != null)
			{
				boolean containingReqSeq = false;
				for(Pair<Integer,String> pair : reqStrList)
				{
					if(token[pair.getFirst()].equalsIgnoreCase(pair.getSecond()))
						containingReqSeq = true;
				}
				if(containingReqSeq == false)
					continue;
			}
			
			float score = Float.parseFloat(token[scoreCol]); 
			if(pepCol >= 0)
			{
				String pep = token[pepCol];
				
				// if there are flanking amino acids (e.g. R.ACDEFK.G), remove them
				int firstDotIndex = pep.indexOf('.');
				int lastDotIndex = pep.lastIndexOf('.');
				if(firstDotIndex < lastDotIndex)
					pep = pep.substring(firstDotIndex+1, lastDotIndex);
				pep = pep.toUpperCase();

//				int ntt = 0;
//				if(file.getName().contains("Tryp"))
//				{
//					if(token[pepCol].charAt(0) == '.' || token[pepCol].charAt(0) == 'K' || token[pepCol].charAt(0) == 'R')
//						ntt++;
//					if(token[pepCol].charAt(lastDotIndex-1) == '.' || token[pepCol].charAt(lastDotIndex-1) == 'K' || token[pepCol].charAt(lastDotIndex-1) == 'R')
//						ntt++;
//					if(ntt == 0)
//						continue;
//				}
//				else if(file.getName().contains("LysN"))
//				{
//					if(token[pepCol].charAt(firstDotIndex+1) == 'K')
//						ntt++;
//					if(token[pepCol].charAt(token[pepCol].length()-1) == '.' || token[pepCol].charAt(token[pepCol].length()-1) == 'K')
//						ntt++;
//					if(ntt == 0)
//						continue;
//				}
				
				
				Float prevScore = prevScoreTable.get(pep);
				if(prevScore == null || (isGreaterBetter && score > prevScore) || (!isGreaterBetter && score < prevScore))
					prevScoreTable.put(pep, score);
			}
			else
				scoreList.add(score);
		}
		if(pepCol >= 0)
		{
			Iterator<Entry<String, Float>> itr = prevScoreTable.entrySet().iterator();
			while(itr.hasNext())
			{
				Entry<String, Float> entry = itr.next();
				list.add(new ScoredString(entry.getKey(), entry.getValue()));
			}
		}
		else
		{
			for(Float score : scoreList)
				list.add(new ScoredString(null, score));
		}
		return list;
	}
	
	public static TreeMap<Float,Float> printROCCurve(String identifier, ArrayList<Float> target, ArrayList<Float> decoy, 
			boolean isGreaterBetter, boolean isConcatenated)
	{
		TreeMap<Float,Float> fdrMap = new TreeMap<Float,Float>();
		float[] fdrInterested = {0.01f};
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
		
		ArrayList<Pair<Float,Integer>> curve = new ArrayList<Pair<Float,Integer>>();
		ArrayList<Float> threshold = new ArrayList<Float>();
		
		float prevFDR = -1;
		int prevNumber = -1;
		int targetIndex = 0;
		
		float[] fdrToBeReported = new float[fdrInterested.length];
		float[] thresholds = new float[fdrInterested.length];
		int[] numTargetMatches = new int[fdrInterested.length];
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
				if(!isConcatenated)
					fdr = decoyIndex/(float)targetIndex;
				else
					fdr = (2*decoyIndex)/(float)targetIndex;

				fdrMap.put(decoyScore, fdr);
				
				if(fdr > prevFDR && targetIndex > prevNumber && fdr <= FDR_REPORT_THRESHOLD)
				{
					curve.add(new Pair<Float,Integer>(fdr, targetIndex));
					threshold.add(decoyScore);
					prevFDR = fdr;
					prevNumber = targetIndex;
				}
				for(int i=0; i<fdrInterested.length; i++)
				{
					if(fdr > fdrInterested[i])
						continue;
					if(targetIndex > numTargetMatches[i])
					{
						fdrToBeReported[i] = fdr;
						thresholds[i] = decoyScore;
						numTargetMatches[i] = targetIndex;
					}
				}
			}
		}
		
		
		if(curve.size() > 0)
		{
			System.out.print(identifier+"_FDR=["+curve.get(0).getFirst());
			for(int i=1; i<curve.size(); i++)
				System.out.print(","+curve.get(i).getFirst());
			System.out.println("];");
			System.out.print(identifier+"_NUM=["+curve.get(0).getSecond());
			for(int i=1; i<curve.size(); i++)
				System.out.print(","+curve.get(i).getSecond());
			System.out.println("];");
			System.out.print(identifier+"_THR=["+threshold.get(0));
			for(int i=1; i<threshold.size(); i++)
				System.out.print(","+threshold.get(i));
			System.out.println("];");
			for(int i=0; i<fdrInterested.length; i++)
			{
				System.out.println("%"+fdrToBeReported[i]+" "+thresholds[i]+" "+numTargetMatches[i]);
			}
			
		}
		
		return fdrMap;
	}
}
