package edu.ucsd.msjava.misc;

import java.io.*;
import java.math.*;
import java.util.*;

import org.systemsbiology.jrap.stax.*;

import edu.ucsd.msjava.msdbsearch.*;
import edu.ucsd.msjava.msgf.*;
import edu.ucsd.msjava.msscorer.*;
import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.msutil.Modification.Location;
import edu.ucsd.msjava.parser.*;


public class Chores {
	public static void main(String argv[]) throws Exception
	{
//		System.out.println(Enzyme.TRYPSIN.getProbCleavageSites());
//		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
//		DBScanner.setAminoAcidProbabilities("/home/sangtaekim/Research/Data/IPI/IPI_human_3.79.fasta", aaSet);
//		System.out.println(Enzyme.TRYPSIN.getProbCleavageSites());
//		System.out.format("%.2f", 10.346f);
//		System.out.println(new IonType.PrefixIon(1,19));
//		long time = System.currentTimeMillis();
//		float value = 1;
//		for(int i=0; i<1000000000; i++)
//			for(int j=0; j<1000000000; j++)
//				value += Math.log(2.3f);
//		System.out.println(System.currentTimeMillis()-time);
//		System.out.println(IonType.X.getOffset());
//		System.out.println(Character.isUpperCase(']'));
//		System.out.println(Float.MIN_NORMAL);
//		System.out.println(Double.MIN_NORMAL);
//		System.out.println(Enzyme.LysN.getPeptideCleavageEfficiency());
//		System.out.println(Enzyme.LysN.getNeighboringAACleavageEffiency());
//		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
//		aaSet.registerEnzyme(Enzyme.LysC);
//		System.out.println(aaSet.getNeighboringAACleavageCredit());
//		System.out.println(aaSet.getNeighboringAACleavagePenalty());
//		System.out.println(aaSet.getPeptideCleavageCredit());
//		System.out.println(aaSet.getPeptideCleavagePenalty());
//		efdrTest();
//		System.out.println((Composition.N15-Composition.N)*2+(Composition.C13-Composition.C)*6);
//		System.out.println((Composition.N15-Composition.N)*4+(Composition.C13-Composition.C)*6);
//		System.out.println(Composition.OFFSET_Y);
//		System.out.println(Runtime.getRuntime().availableProcessors());
//		System.out.println(Composition.getMass("H-2O-1"));
//		printMaxScanNum();
//		aaSet.registerEnzyme(Enzyme.TRYPSIN);
//		System.out.println(aaSet.getNeighboringAACleavageCredit()+" "+aaSet.getNeighboringAACleavagePenalty());
//		System.out.println(aaSet.getPeptideCleavageCredit()+" "+aaSet.getPeptideCleavagePenalty());
//		testSpecKey();
//		String pep = "-.asdf.AA";
//		System.out.println(pep.matches(".\\..+\\.."));
//		System.out.println((byte)200);
//		simpleTest();
//		System.out.println("1");
//		System.out.println(Composition.getMass("H-2O-1"));
//		System.out.println(Composition.getMass("C2H2"));
//		Modification phospho = Modification.get("Phosphorylation");
//		System.out.println(phospho.getAccurateMass());
//		recombCPReg();
		
//		ArrayList<Modification.Instance> mods = new ArrayList<Modification.Instance>();
//		mods.add(new Modification.Instance(Modification.get("Carbamidomethylation"), 'C').fixedModification());
//		mods.add(new Modification.Instance(Modification.get("PyroCarbamidomethyl"), 'C', Location.N_Term));
//		mods.add(new Modification.Instance(Modification.get("Oxidation"), 'M', Location.Anywhere));
//		mods.add(new Modification.Instance(Modification.get("Acetylation"), '*', Location.N_Term));
//		mods.add(new Modification.Instance(Modification.get("PyrogluQ"), 'Q', Location.N_Term));
//		mods.add(new Modification.Instance(Modification.get("PyrogluE"), 'E', Location.N_Term));
//		
//		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet(mods);
//		aaSet.printAASet();

//		String s = "SEQ=(Q,-17.026549)YWEYQFQHQPSQEE(C,57.02146)EGSSLSAVFEHFAMMQR";
//		String[] token = s.split("[(,)]");
//		for(int i=0; i<token.length; i++)
//			System.out.println(i+"\t"+token[i]);
		
//		String s = "SEQ=(A,42.0106)AAAA(C,39.99)AVGPGAGGAGSAVPGGAGP(C,57.02146)(C,57.02146)ATVSVFPGAR";
//		s = s.replaceAll("\\(C,39\\.99\\)", "\\(C,-17\\.026549\\)");
//		s = s.replaceAll("\\(C,57\\.02146\\)", "C");
//		System.out.println(s);
		
//		String buf = "SCANS=953-939";
//		System.out.println(buf.matches(".+=\\d+-\\d+"));
//		int startScanNum = Integer.parseInt(buf.substring(buf.indexOf('=')+1, buf.lastIndexOf('-')));
//		int endScanNum = Integer.parseInt(buf.substring(buf.lastIndexOf('-')+1));
//		System.out.println(startScanNum + " " + endScanNum);

//		System.out.println(Composition.C13 - Composition.C);
//		System.out.println(Composition.N15 - Composition.N);
//		System.out.println(System.getProperty("user.home"));
		
		String compStr = "C2H3NOFe";
		System.out.println(compStr + ": " + Composition.getMass(compStr));
	}

	public static void recombCPReg() throws Exception
	{
		File regList = new File("/Users/sangtaekim/Dropbox/Documents/RECOMB-CP 2012/regList.txt");
		File submissions = new File("/Users/sangtaekim/Dropbox/Documents/RECOMB-CP 2012/submissions.txt");
		
		HashSet<String> names = new HashSet<String>();
		
		String s;
		BufferedLineReader in = new BufferedLineReader(regList.getPath());
		while((s=in.readLine()) != null)
		{
			names.add(s.toLowerCase());
		}
		in.close();
					
		in = new BufferedLineReader(submissions.getPath());
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length != 4)
			{
				System.err.println("Error: " + s);
				System.exit(-1);
			}
			String authors = token[1];
			String type = token[3];
			if(!type.equalsIgnoreCase("Abstract"))
				continue;
			String[] aToken = authors.split("\\s+");
			boolean isRegistered = false;
			for(String name : aToken)
			{
				for(String r : names)
				{
					if(r.contains(name.toLowerCase()))
					{
						System.out.println("***" + name + ":" + r);
						isRegistered = true;
					}
				}
			}
			if(isRegistered)
				System.out.println(s);
		}
		
	}
	
	public static void simpleTest() throws Exception
	{
		int specSize = 728178;
		int fromIndexGlobal = 0;
		int numThreads = 8;
		int numSpecScannedTogether = 116510;
		
		while(true)
		{
			if(fromIndexGlobal >= specSize)
				break;
			int toIndexGlobal = Math.min(specSize, fromIndexGlobal+numSpecScannedTogether);
			System.out.println("=========== " + fromIndexGlobal + "," + toIndexGlobal + " =============");
			int size = toIndexGlobal - fromIndexGlobal;
			int subListSize = size/numThreads;
			int residue = size % numThreads;
			
			int[] startIndex = new int[numThreads];
			int[] endIndex = new int[numThreads];
			
			for(int i=0; i<numThreads; i++)
			{
				startIndex[i] =  i > 0 ? endIndex[i-1] : fromIndexGlobal;
				endIndex[i] = startIndex[i] + subListSize + (i < residue ? 1 : 0);
				
				System.out.println(startIndex[i]+","+endIndex[i]);
			}
			fromIndexGlobal += numSpecScannedTogether;
		}
		
	}
	
	public static void printMaxScanNum() throws Exception
	{
		String mzXMLFileName = "/home/sangtaekim/Test/CCMS/3297b97db35241ba8547906b22377869/spectrum/00000.mzXML";
		MzXMLSpectraMap map = new MzXMLSpectraMap(mzXMLFileName);
		System.out.println(map.getMaxScanNumber());
	}
	
	public static void compareFiles() throws Exception
	{
		String mzXMLFileName = "/home/sangtaekim/Test/CCMS/3297b97db35241ba8547906b22377869/spectrum/00000.mzXML";
		String mgfFileName = "/home/sangtaekim/Test/CCMS/3297b97db35241ba8547906b22377869/spectrum/00000.mgf";
		
		SpectraIterator itr = new SpectraIterator(mgfFileName, new MgfSpectrumParser());
		HashSet<Integer> mgfScanNum = new HashSet<Integer>();
		while(itr.hasNext())
		{
			mgfScanNum.add(itr.next().getScanNum());
		}
		
		MzXMLSpectraIterator itr2 = new MzXMLSpectraIterator(mzXMLFileName);
		HashSet<Integer> mzXMLScanNum = new HashSet<Integer>();
		while(itr2.hasNext())
		{
			mzXMLScanNum.add(itr2.next().getScanNum());
		}
		
		int mzXMLOnly = 0;
		for(int scanNum : mzXMLScanNum)
		{
			if(!mgfScanNum.contains(scanNum))
			{
				System.out.println("MzXMLOnly: " + scanNum);
				mzXMLOnly++;
			}
		}

		int mgfOnly = 0;
		for(int scanNum : mgfScanNum)
		{
			if(!mzXMLScanNum.contains(scanNum))
			{
				System.out.println("MgfOnly: " + scanNum);
				mgfOnly++;
			}
		}
		
		System.out.println("MzXMLOnly: " + mzXMLOnly);
		System.out.println("MgfOnly: " + mgfOnly);
		
	}
	public static void efdrTest() throws Exception
	{
		double[] specProb = new double[136964];
		Random rand = new Random();
		for(int i=0; i<specProb.length; i++)
		{
			specProb[i] = rand.nextDouble();
//			Arrays.sort(specProb[i]);
		}
		
		long time = System.currentTimeMillis();
		Arrays.sort(specProb);
		System.out.println("Sorting: " + (System.currentTimeMillis()-time));
		time = System.currentTimeMillis();
		double cumP = 0;
		for(int i=0; i<specProb.length; i++)
		{
			double probCorr = 1.-specProb[i];
			double pValue = MSGFDBResultGenerator.DBMatch.getPValue(specProb[i], 3000000);
			cumP += pValue;
			double eTD = i+1-cumP;
			double eDD = cumP;
			for(int j=1; j<specProb.length; j++)
			{
				eDD += specProb[j];
			}
			double eFDR = eDD/eTD;
			double dummy = eFDR;
		}
		System.out.println("Time: " + (System.currentTimeMillis()-time));
		
	}
	
	public static void combination(int n, int r) throws Exception
	{
		for(LinkedList<Integer> ins : getCombinations(n,r))
			System.out.println(ins);
	}

	public static ArrayList<LinkedList<Integer>> getCombinations(int n, int r) throws Exception
	{
		ArrayList<LinkedList<Integer>> nHr = new ArrayList<LinkedList<Integer>>();
		
		if(r==1)
		{
			for(int i=0; i<n; i++)
			{
				LinkedList<Integer> newIns = new LinkedList<Integer>();
				newIns.add(i);
				nHr.add(newIns);
			}
		}
		else if(r>1)
		{
			ArrayList<LinkedList<Integer>> nHrMinus1 = getCombinations(n, r-1);
			for(LinkedList<Integer> ins : nHrMinus1)
			{
				int prevLargest = ins.getLast();
				for(int i=prevLargest; i<n; i++)
				{
					LinkedList<Integer> newIns = new LinkedList<Integer>(ins);
					newIns.add(i);
					nHr.add(newIns);
				}
			}
		}
		
		return nHr;
	}
	
	public static void mzXMLLoadingTest() throws Exception
	{
		long time = System.currentTimeMillis();
		String fileName = "/home/sangtaekim/Research/Data/HeckWhole/Spectra/090121_NM_Trypsin_20.mzXML";
		MzXMLSpectraMap map = new MzXMLSpectraMap(fileName);
		System.out.println("LoadingTime: " + (System.currentTimeMillis()-time));
		map.getSpectrumBySpecIndex(100);
		System.out.println("Time: " + (System.currentTimeMillis()-time));
	}
	
	public static void extractRECOMBCPEmails() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Documents/RECOMB-CP/participantsEmails.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.matches(".+@.+\\..+"))
				System.out.println(s.split("\\s+")[0]);
		}
	}
	public static void extractRECOMBBEAbstracts() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Desktop/abstracts.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		int lineNum=0;
		Hashtable<String, String> entries = new Hashtable<String,String>();
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("ID"))
				continue;
			lineNum++;
			String[] token = s.split("\t");
			if(token.length < 12)
				continue;
			String email = token[7];
			if(!email.contains("sangtae") && !email.contains("chorta"))
				entries.put(email, s);
		}
		
		ArrayList<String> list = new ArrayList<String>();
		for(String e : entries.values())
			list.add(e);
		Collections.sort(list);
		for(String e : list)
		{
			String[] token = e.split("\t");
			System.out.println("****************************");
			System.out.println("Date Submitted: " + token[1]);
			System.out.println("Category: " + token[2]);
			System.out.println("Title: " + token[3]);
			System.out.println("Author: " + token[4] + " " + token[5]);
			System.out.println("Co-Authors: " + token[6]);
			System.out.println("Email: " + token[7]);
			System.out.println("Affiliation: " + token[8]);
			System.out.println("Topic: " + token[9]);
			System.out.println("Keywords: " + token[10]);
			System.out.println();
			System.out.println("Abstract:");
			for(int i=11; i<token.length; i++)
				System.out.print(token[i]+"\t");
			System.out.println();
			System.out.println();
		}
	}
	
	public static void extractEmails() throws Exception
	{
		String fileName = "/Users/sangtaekim/Documents/RECOMB-CP/asms09.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\\s+");
			for(String t : token)
			{
				if(t.contains("@") && t.contains("."))
				{
					System.out.println(t);
				}
			}
		}
	}
}
