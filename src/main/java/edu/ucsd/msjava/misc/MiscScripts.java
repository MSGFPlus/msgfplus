package edu.ucsd.msjava.misc;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.Map.Entry;


import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.ScoredString;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.parser.MgfSpectrumParser;
import edu.ucsd.msjava.suffixarray.SuffixArray;
import edu.ucsd.msjava.suffixarray.SuffixArraySequence;

public class MiscScripts {
	public static void main(String argv[]) throws Exception
	{
//		makeMatchedProteinDB();
//		processShortAgilentPeptides();
//		suffixArrayTest();
//		processShortAgilentPeptidesAll(5);
//		countSpectra(0, 10000, 2, "/home/sangtaekim/Research/Data/Heck/annotatedHECK_CID_LysN.mgf");
		peptideTest();
	}

	public static void peptideTest() 
	{
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		Peptide pep = aaSet.getPeptide("GCDEFG");
		if(pep != null)
			for(AminoAcid aa : pep)
				System.out.println(aa.getResidueStr()+"\t"+aa.getMass());
		else
			System.out.println("NULL");
		
	}
	public static void processShortAgilentPeptidesAll(int length) throws Exception
	{
		System.out.print("Length\tMSGFThreshold\t#Positive\t#TargetMatch\tRatioTargetMatch");
		System.out.println("\t#DecoyMatch\tFDR\t#ExprDBMatch\t#ExprDecoyMatch\tExprFDR\t" +
				"#TargetPep\t#DecoyPep\tFDRPep\t#ExprPep\t#ExprDecoyPep\tFDRExprPep" +
				"\tNumDeNovoPeptides");

		String deNovoScoreFileName = System.getProperty("user.home")+"/Research/Data/SignalPeptides/agilentDeNovoScores.txt";
		BufferedLineReader in2 = new BufferedLineReader(deNovoScoreFileName);
		String s;
		in2.readLine();
		int[] scoreHist = new int[100];
		int[] numRecs = new int[100];
		while((s=in2.readLine()) != null)
		{
			String[] token = s.split("\t");
			assert(token.length == 4);
			int l = Integer.parseInt(token[1]);
			int score = Integer.parseInt(token[2]);
			int numDeNovo = Integer.parseInt(token[3]);
			if(l == length || l == '*')
			{
				scoreHist[score]++;
				numRecs[score]+=numDeNovo;
			}
		}
		in2.close();
		
		String fileName = System.getProperty("user.home")+"/Research/Data/SignalPeptides/agilentDeNovoAll.txt";
		
		ArrayList<ScoredString> msgfTarget = new ArrayList<ScoredString>();
		ArrayList<ScoredString> msgfDecoy = new ArrayList<ScoredString>();
		ArrayList<ScoredString> msgfExprTarget = new ArrayList<ScoredString>();
		ArrayList<ScoredString> msgfExprDecoy = new ArrayList<ScoredString>();

		BufferedLineReader in = new BufferedLineReader(fileName);
		int datasetIndex = -1;
		ArrayList<ScoredString> curList = null;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#"))
			{
				datasetIndex++;
				if(datasetIndex == 0)
					curList = msgfTarget;
				else if(datasetIndex == 1)
					curList = msgfDecoy;
				else if(datasetIndex == 2)
					curList = msgfExprTarget;
				else if(datasetIndex == 3)
					curList = msgfExprDecoy;
				in.readLine();
				in.readLine();
			}
			String[] token = s.split("\t");
			assert(token.length == 9);
			int scanNum = Integer.parseInt(token[1]);
			String peptide = token[2].substring(token[2].indexOf('.')+1, token[2].lastIndexOf('.'));
			int msgfScore = Integer.parseInt(token[6]);
			float specProb = Float.parseFloat(token[8]);
			curList.add(new ScoredString(scanNum+":"+peptide+":"+specProb, msgfScore));
		}
		
		for(int threshold=0; threshold<=42; threshold++)
		{
			int numPositives = 0;
			float avgNumDeNovoPeptides=0;
			for(int t = threshold; t<scoreHist.length; t++)
			{
				numPositives += scoreHist[t];
				avgNumDeNovoPeptides += numRecs[t];
			}
			avgNumDeNovoPeptides /= numPositives;
			
			HashSet<String> scanNumSet = new HashSet<String>();
			HashSet<String> pepSet = new HashSet<String>();
			for(ScoredString ss : msgfTarget)
			{
				String[] token = ss.getStr().split(":");
				String scanNum = token[0];
				String peptide = token[1];
				if(peptide.length() == length && ss.getScore() >= threshold)
				{
					scanNumSet.add(scanNum);
					pepSet.add(peptide);
				}
			}
			int numTargetSpectra = scanNumSet.size();
			int numTargetPep = pepSet.size();
			
			scanNumSet = new HashSet<String>();
			pepSet = new HashSet<String>();
			for(ScoredString ss : msgfDecoy)
			{
				String[] token = ss.getStr().split(":");
				String scanNum = token[0];
				String peptide = token[1];
				if(peptide.length() == length && ss.getScore() >= threshold)
				{
					scanNumSet.add(scanNum);
					pepSet.add(peptide);
				}
			}
			int numDecoySpectra = scanNumSet.size();
			int numDecoyPep = pepSet.size();
			
			scanNumSet = new HashSet<String>();
			pepSet = new HashSet<String>();
			for(ScoredString ss : msgfExprTarget)
			{
				String[] token = ss.getStr().split(":");
				String scanNum = token[0];
				String peptide = token[1];
				if(peptide.length() == length && ss.getScore() >= threshold)
				{
					scanNumSet.add(scanNum);
					pepSet.add(peptide);
				}
			}
			int numExprTargetSpectra = scanNumSet.size();
			int numExprTargetPep = pepSet.size();
	
			scanNumSet = new HashSet<String>();
			pepSet = new HashSet<String>();
			for(ScoredString ss : msgfExprTarget)
			{
				String[] token = ss.getStr().split(":");
				String scanNum = token[0];
				String peptide = token[1];
				if(peptide.length() == length && ss.getScore() >= threshold)
				{
					scanNumSet.add(scanNum);
					pepSet.add(peptide);
				}
			}
			int numExprDecoySpectra = scanNumSet.size();
			int numExprDecoyPep = pepSet.size();
			
			System.out.print(length+"\t"+threshold+"\t"+numPositives+"\t"+numTargetSpectra+"\t"+(numTargetSpectra/(float)numPositives));
			System.out.print("\t"+numDecoySpectra+"\t"+(numDecoySpectra/(float)numTargetSpectra));
			System.out.print("\t"+numExprTargetSpectra+"\t"+numExprDecoySpectra+"\t"+(numExprDecoySpectra/(float)numExprTargetSpectra));
			System.out.print("\t"+numTargetPep+"\t"+numDecoyPep+"\t"+(numDecoyPep/(float)numTargetPep));
			System.out.print("\t"+numExprTargetPep+"\t"+numExprDecoyPep+"\t"+(numExprDecoyPep/(float)numExprTargetPep));
			System.out.println(avgNumDeNovoPeptides);
		}
	}
	
	public static void countSpectra(float minMass, float maxMass, int charge, String fileName) throws Exception
	{
//		String fileName = "/Users/sangtaekim/Research/Data/AgilentQTOF/agilentQTOFAll.mgf";
		SpectraIterator iterator = new SpectraIterator(fileName, new MgfSpectrumParser());
		int numSpecs = 0;
		HashSet<String> pepSet = new HashSet<String>();
		while(iterator.hasNext())
		{
			Spectrum spec = iterator.next();
			if(spec.getCharge() != charge && spec.getCharge() != 0)
				continue;
			spec.getPrecursorPeak().setCharge(2);
			float parentMass = spec.getParentMass();
			if(parentMass < minMass || parentMass > maxMass)
				continue;
			numSpecs++;
			if(spec.getAnnotationStr() != null)
				pepSet.add(spec.getAnnotationStr());
		}
		System.out.println("NumSpecs: " + numSpecs);
		System.out.println("NumPeptides: " + pepSet.size());
	}
	

	public static void processShortAgilentPeptides() throws Exception
	{
		System.out.print("Length\tMSGFThreshold\tDeltaThreshold\t#Positive\t#TargetMatch\tRatioTargetMatch");
		System.out.println("\t#DecoyMatch\tFDR\t#ExprDBMatch\t#ExprDecoyMatch\tExprFDR\t" +
				"#TargetPep\t#DecoyPep\tFDRPep\t#ExprPep\t#ExprDecoyPep\tFDRExprPep" +
				"\tNumDeNovoPeptides");
		String targetFileName = System.getProperty("user.home")+"/Research/Data/SignalPeptides/agilentDeNovoFwd.txt";
		String decoyFileName = System.getProperty("user.home")+"/Research/Data/SignalPeptides/agilentDeNovoRev.txt";
		int length = 6;
//		int deltaThreshold = 0;
		int msgfThreshold = 30;
//		for(int msgfThreshold = 0; msgfThreshold <= 43; msgfThreshold++)
		for(int deltaThreshold = 0; deltaThreshold <= 10; deltaThreshold++)
		{
			String target = processShortAgilentPeptides(targetFileName, length, msgfThreshold, deltaThreshold);
			String decoy = processShortAgilentPeptides(decoyFileName, length, msgfThreshold, deltaThreshold);
			String numDecoyMatches = decoy.split("\t")[4];
			int numDecoyPep = Integer.parseInt(decoy.split("\t")[6]);
			int numExprDecoyMatch = processShortAgilentDecoyExpr(length, msgfThreshold, false);
			if(numExprDecoyMatch == 0)
				continue;
			String[] token = target.split("\t");
			for(int i=0; i<5; i++)
				System.out.print(token[i]+"\t");
			
			System.out.print(Integer.parseInt(token[4])/(float)Integer.parseInt(token[3])+"\t");
			System.out.print(numDecoyMatches+
					"\t" + Integer.parseInt(numDecoyMatches)/(float)Integer.parseInt(token[4]));
			System.out.print("\t"+ token[5] + "\t" + numExprDecoyMatch +
					"\t"+numExprDecoyMatch/(float)Integer.parseInt(token[5]));

			int numTargetPep = Integer.parseInt(token[6]);
			int numExprPep = Integer.parseInt(token[7]);
			int numExprDecoyPep = processShortAgilentDecoyExpr(length, msgfThreshold, true);
			float numDeNovo = Float.parseFloat(token[8]);
			System.out.println("\t"+numTargetPep+"\t"+numDecoyPep+"\t"+numDecoyPep/(float)numTargetPep+"\t"+
					numExprPep+"\t"+numExprDecoyPep+"\t"+numExprDecoyPep/(float)numExprPep+"\t"+numDeNovo);
		}
	}
	
	private static int processShortAgilentDecoyExpr(int length, int threshold, boolean uniquePep) throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/SignalPeptides/agilentDeNovoExprDecoy.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		in.readLine();	// annotation
		//int count = 0;
		HashSet<String> titleSet = new HashSet<String>();
		HashSet<String> pepSet = new HashSet<String>();
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length != 9)
				continue;
			int msgfScore = Integer.parseInt(token[6]);
			if(msgfScore < threshold)
				continue;
			String pep = token[2].substring(token[2].indexOf('.')+1, token[2].lastIndexOf('.'));
			if(pep.length() != length)
				continue;
			titleSet.add(token[0]);
			pepSet.add(pep);
		}
		if(uniquePep)
			return pepSet.size();
		else
			return titleSet.size();
	}
	
	private static String processShortAgilentPeptides(String fileName, int length, int msgfThreshold, int deltaThreshold) throws Exception
	{
		String s;
		BufferedLineReader lineReader = new BufferedLineReader(fileName);
		int numPositive = 0;
		int numTarget = 0;
		int numExpr = 0;
		int sumRecs = 0;
		int numRecs = 0;
		boolean passThreshold = false;
		boolean lengthMatch = false;
		boolean targetMatch = false;
		boolean exprMatch = false;
		String pep = null;
		String exprPep = null;
		HashSet<String> pepSet = new HashSet<String>();
		HashSet<String> exprPepSet = new HashSet<String>();
		while((s = lineReader.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length >= 5)
			{
				if(passThreshold && lengthMatch)
				{
					numPositive++;
					if(targetMatch)
						numTarget++;
					if(exprMatch)
						numExpr++;
					if(pep != null)
						pepSet.add(pep);
					if(exprPep != null)
						exprPepSet.add(exprPep);
					sumRecs += numRecs;
				}
				passThreshold = false;
				lengthMatch = false;
				targetMatch = false;
				exprMatch = false;
				pep = null;
				exprPep = null;
				numRecs = 0;
					
				//String title = token[0];
				//float precurMass = Float.parseFloat(token[1]);
				//int charge = Integer.parseInt(token[2]);
				int msgfScore = Integer.parseInt(token[3]);
				int deltaScore = Integer.parseInt(token[4]);
				if(msgfScore >= msgfThreshold && deltaScore >= deltaThreshold)
				{
					passThreshold = true;
				}
			}
			else if(token.length == 3)
			{
				if(token[0].length() == length && passThreshold)
				{
					numRecs++;
					lengthMatch = true;
					if(token[1].equalsIgnoreCase("1"))
					{
						targetMatch = true;
						pep = token[0];
					}
					if(token[2].equalsIgnoreCase("1"))
					{
						exprMatch = true;
						exprPep = token[0];
					}
				}
			}
		}
		if(passThreshold && lengthMatch)
		{
			numPositive++;
			if(targetMatch)
				numTarget++;
			if(exprMatch)
				numExpr++;
			if(pep != null)
				pepSet.add(pep);
			if(exprPep != null)
				exprPepSet.add(exprPep);
			sumRecs += numRecs;
		}

		return length+"\t"+msgfThreshold+"\t"+deltaThreshold+"\t"+numPositive+"\t"+numTarget+
		"\t"+numExpr+"\t"+pepSet.size()+"\t"+exprPepSet.size()+"\t"+sumRecs/(float)numPositive;
	}
	
	public static void suffixArrayTest() throws Exception
	{
		int splitNum = 0;
		String dbFileName = System.getProperty("user.home")+"/Research/Data/HumanGenome/translated/HSRM.NCBI36.54.translation."+splitNum+".fasta";
		SuffixArray sa = new SuffixArray(new SuffixArraySequence(dbFileName, edu.ucsd.msjava.sequences.Constants.AMINO_ACIDS_19));
		System.out.println(sa.search("STHQHSENEEFPKL"));
	}
	
	public static void makeMatchedProteinDB() throws Exception
	{
		String dbFileName = System.getProperty("user.home")+"/Research/Data/EColiDB/Ecol_protein_formatted.fasta";
		SuffixArraySequence adapter = new SuffixArraySequence(dbFileName, edu.ucsd.msjava.sequences.Constants.AMINO_ACIDS_19);
		SuffixArray sa = new SuffixArray(adapter);
		
		String specFileName = System.getProperty("user.home")+"/Research/Data/AgilentQTOF/annotatedAgilentQTOF.mgf";
		SpectraIterator iterator = new SpectraIterator(specFileName, new MgfSpectrumParser());

		TreeMap<String, String> proteinDB = new TreeMap<String, String>();
		
		int specNum = 0;
		while(iterator.hasNext())
		{
			Spectrum spec = iterator.next();
			specNum++;
			String peptide = spec.getAnnotationStr();
			if(peptide.length() >= 7)
			{
				ArrayList<String> annotations = sa.getAllMatchingAnnotations(peptide);
				ArrayList<String> proteins = sa.getAllMatchingEntries(peptide);
				assert(annotations.size() == proteins.size());
				for(int i=0; i<annotations.size(); i++)
					proteinDB.put(annotations.get(i), proteins.get(i));
			}
		}
		for(Entry<String, String> entry : proteinDB.entrySet())
		{
			System.out.println(">"+entry.getKey());
			System.out.println(entry.getValue());
		}
	}
}
