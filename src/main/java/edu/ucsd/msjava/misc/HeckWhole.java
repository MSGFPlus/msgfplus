package edu.ucsd.msjava.misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;


import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.IonType;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.SpectraContainer;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.SpectraMap;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.SpectrumAccessorBySpecIndex;
import edu.ucsd.msjava.msutil.SpectrumAnnotator;
import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.parser.MgfSpectrumParser;
import edu.ucsd.msjava.parser.MzXMLSpectraIterator;
import edu.ucsd.msjava.parser.MzXMLSpectraMap;
import edu.ucsd.msjava.parser.MzXMLToMgfConverter;

public class HeckWhole {
	public static void main(String argv[]) throws Exception
	{
//		convertMzXMLIntoMgf();
//		extractGoodSpectra();
//		precursorTest();
//		rankOnePeakTest();
//		checkGoodSpectra("Tryp");
//		extractGoodSpectraFromSum();
//		vennDiagram();
//		test();
//		vennDiagramMascotMSGF();
//		deNovoPerformance();
//		splitMascotResults();
		countEnzymeCleavage("ETD","Tryp");
	}

	public static void splitMascotResults() throws Exception
	{
		String dir = System.getProperty("user.home")+"/Research/Data/HeckWhole/ResultsShabaz";
		String[] methods = {"ETD", "CID"};
		String[] enzymes = {"Tryp", "LysN"};
		for(String enzyme : enzymes)
		{
			for(String method : methods)
			{
				BufferedLineReader in = new BufferedLineReader(dir+"/"+enzyme+"_"+method+".txt");
				PrintStream outTarget = new PrintStream(new BufferedOutputStream(new FileOutputStream(dir+"/"+enzyme+"_"+method+"_GygiTarget.txt")));
				PrintStream outDecoy = new PrintStream(new BufferedOutputStream(new FileOutputStream(dir+"/"+enzyme+"_"+method+"_GygiDecoy.txt")));
				String s;
				String prevTitle = "";
				while((s=in.readLine()) != null)
				{
					if(s.startsWith("#"))
						continue;
					String[] token = s.split("\t");
					assert(token.length == 5);
					String title = token[0];
					if(title.equalsIgnoreCase(prevTitle))
						continue;
					else
						prevTitle = title;
					String protein = token[4];
					PrintStream out;
					if(protein.contains("IPI:REV"))
						out = outDecoy;
					else
						out = outTarget;
					out.println(s);
				}
				outTarget.close();
				outDecoy.close();
			}
		}
	}
	
	public static void precursorOFF() throws Exception
	{
		String[] methods = {"ETD", "CID"};
		String[] enzymes = {"Tryp", "LysN"};
		int[] charges = {2,3};
		int[] chargeReduces = {0, 1};
		
		for(String method : methods)
			for(String enzyme : enzymes)
				for(int charge : charges)
					for(int chargeReduce : chargeReduces)
						precursorOFF(method, enzyme, charge, chargeReduce);
	}
	
	public static void precursorOFF(String method, String enzyme, int charge, int chargeReduce) throws Exception
	{
		if(charge-chargeReduce < 2)
			return;
		String fileName = System.getProperty("user.home")+"/Research/Data/HeckWhole/AnnotatedSpectra/annotatedHeckWhole_"+method+"_"+enzyme+".mgf";
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		int[] hist = new int[401];
		int numSpectra = 0;
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			if(spec.getCharge() != charge)
				continue;
			numSpectra++;
			spec.setRanksOfPeaks();
			float precursorMz = spec.getPrecursorPeak().getMz();
			float neutralParentMass = (precursorMz-(float)Composition.NEUTRON)*spec.getCharge();
			float mass = (neutralParentMass+(charge-chargeReduce)*(float)Composition.NEUTRON)/(charge-chargeReduce);
			for(Peak p : spec)
			{
				if(p.getRank() <= 20)
				{
					float offset = p.getMz() - mass;
					if(offset > -100 && offset < 100)
					{
						int binNum = Math.round(offset/0.5f)+200;
						hist[binNum]++;
					}
				}
			}
		}
		
		System.out.println(method+enzyme);
		System.out.println("#Spectra\t"+numSpectra);
		System.out.println("Charge\t"+charge);
		System.out.println("ConsideredCharge\t"+(charge-chargeReduce));
		System.out.println("Offset\t#Peaks(within rank 20)");
		for(int i=1; i<hist.length-1; i++)
		{
			System.out.format("%.1f\t%f\n", -100+0.5f*i, hist[i]/(float)numSpectra);
		}
		System.out.println();
	}	
	
	private static Hashtable<String, Float[]> mascotThresholds;
	private static Hashtable<String, Float[]> msgfThresholds;

	static {
		Float[]	mascotThresholdCIDTryp = {37.09f, 23.37f, 24.18f};
		Float[]	mascotThresholdETDTryp = {46.69f, 75.14f, 44.2f};
		Float[]	mascotThresholdCIDLysN = {33.93f, 22.35f, 23.01f};
		Float[]	mascotThresholdETDLysN = {35.8f, 30.68f, 22.71f};
		mascotThresholds = new Hashtable<String, Float[]>();		
		mascotThresholds.put("CIDTryp", mascotThresholdCIDTryp);
		mascotThresholds.put("ETDTryp", mascotThresholdETDTryp);
		mascotThresholds.put("CIDLysN", mascotThresholdCIDLysN);
		mascotThresholds.put("ETDLysN", mascotThresholdETDLysN);
		
//		// After-training: for Confident, aa freq modified
		Float[]	thresholdCIDTryp = {1.216313E-10f, 8.098124E-11f, 2.1928117E-11f};
		Float[]	thresholdETDTryp = {8.2302345E-11f, 1.6922793E-12f, 3.0688185E-13f};
		Float[]	thresholdSumTryp = {1.2440679E-10f, 4.8907045E-11f, 1.927252E-11f};
		Float[]	thresholdCIDLysN = {8.042938E-11f, 4.9014993E-11f, 9.152313E-13f};
		Float[]	thresholdETDLysN = {5.6160354E-11f, 4.0307476E-11f, 1.5572712E-11f};
		Float[]	thresholdSumLysN = {6.9574575E-11f, 5.872495E-11f, 1.467543E-11f};
		
//		// After-training: for Confident
//		Float[]	thresholdCIDTryp = {5.2116884E-11f, 2.5660874E-11f, 3.7908413E-12f};
//		Float[]	thresholdETDTryp = {3.3350836E-11f, 6.884241E-13f, 4.9390322E-14f};
//		Float[]	thresholdSumTryp = {5.5019527E-11f, 2.4261168E-11f, 2.7325312E-12f};
//		Float[]	thresholdCIDLysN = {3.4036378E-11f, 1.7352194E-11f, 7.364699E-14f};
//		Float[]	thresholdETDLysN = {2.9340887E-11f, 2.6010118E-11f, 8.103846E-12f};
//		Float[]	thresholdSumLysN = {2.6453065E-11f, 3.1964584E-11f, 4.606455E-12f};

		//		// Pre-training
//		Float[]	thresholdCIDTryp = {5.4535914E-11f, 2.6583753E-11f, 1.8636835E-13f};
//		Float[]	thresholdETDTryp = {3.3466424E-11f, 9.946741E-13f, 6.387666E-14f};
//		Float[]	thresholdSumTryp = {5.6320247E-11f, 2.3743533E-11f, 1.8667672E-11f};
//		Float[]	thresholdCIDLysN = {3.1250006E-11f, 1.3972096E-11f, 1.3098458E-12f};
//		Float[]	thresholdETDLysN = {2.8747773E-11f, 3.3369512E-11f, 1.41990265E-11f};
//		Float[]	thresholdSumLysN = {2.4446337E-11f, 2.8512514E-11f, 9.074788E-12f};

		// Old thresholds
//		Float[]	thresholdCIDTryp = {4.9982934E-11f, 2.5958782E-11f, 1.0835081E-11f};
//		Float[]	thresholdETDTryp = {2.3154301E-11f, 7.760064E-13f, 2.299933E-13f};
//		Float[]	thresholdSumTryp = {5.6320247E-11f, 2.3743533E-11f, 1.8667672E-11f};
//		Float[]	thresholdCIDLysN = {3.2987557E-11f, 1.5884185E-11f, 1.615089E-12f};
//		Float[]	thresholdETDLysN = {2.2711777E-11f, 2.1367989E-11f, 9.767584E-12f};
//		Float[]	thresholdSumLysN = {2.4446337E-11f, 2.8512514E-11f, 9.074788E-12f};
		
		msgfThresholds = new Hashtable<String, Float[]>();		
		msgfThresholds.put("CIDTryp", thresholdCIDTryp);
		msgfThresholds.put("ETDTryp", thresholdETDTryp);
		msgfThresholds.put("SumTryp", thresholdSumTryp);
		msgfThresholds.put("CIDLysN", thresholdCIDLysN);
		msgfThresholds.put("ETDLysN", thresholdETDLysN);
		msgfThresholds.put("SumLysN", thresholdSumLysN);
	}
	
	public static void deNovoPerformance() throws Exception
	{
		String[] methods = {/*"CID", "ETD", */"Sum"};
		String[] enzyme = {"Tryp", "LysN"};
		int[] charges = {2,3};
		for(String m : methods)
			for(String e : enzyme)
				for(int c : charges)
					deNovoPerformance(m, e, c);
	}
	
	public static void deNovoPerformance(String method, String enzyme, int charge) throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/HeckWhole/AnnotatedSpectra/annotatedHeckWhole_"+method+"_"+enzyme+".mgf";
		boolean isSum = false;
		if(method.equals("Sum"))
		{
			isSum = true;
			fileName = System.getProperty("user.home")+"/Research/Data/HeckWhole/AnnotatedSpectra/annotatedHeckWholeSum_CID"+"_"+enzyme+".mgf";
		}
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		int numSpecs = 0;
		int numCorrDeNovo = 0;
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			if(spec.getCharge() != charge)
				continue;
			numSpecs++;
			String title = spec.getTitle();
			String[] token = title.split(":");
			int msgfScore;
			if(!isSum)
				msgfScore = Integer.parseInt(token[3]);
			else
				msgfScore = Integer.parseInt(token[4]);
			int peptideScore;
			if(!isSum)
				peptideScore = Integer.parseInt(token[4]);
			else
				peptideScore = Integer.parseInt(token[5]);
			if(msgfScore == peptideScore)
				numCorrDeNovo++;
		}
		System.out.println((method+enzyme)+"\t"+charge+"\t"+numSpecs+"\t"+numCorrDeNovo+"\t"+(numCorrDeNovo/(float)numSpecs));
	}
	public static void vennDiagramMascotMSGF() throws Exception
	{
//		System.out.println("FragMethod\tEnzyme\tCharge\tMascotOnly\tMSGFOnly\tIntersection\tConflict");
		System.out.println("FileName,ScanNum,Peptide,MascotScore,MSGFPValue");
		String[] enzyme = {"Tryp", "LysN"};
		String[] methods = {"CID","ETD"};
		int[] charges = {2,3};
		for(String e : enzyme)
			for(String m : methods)
				for(int c : charges)
					vennDiagramMascotMSGF(e, m, c);
	}
	
	public static void vennDiagramMascotMSGF(String enzyme, String method, int charge) throws Exception
	{
		String db = "Target";
		
		String mascotFileName = System.getProperty("user.home")+"/Research/Data/HeckWhole/ResultsShabaz/"+enzyme+"_"+method+"_"+db+".txt";
		Hashtable<String,String> mascotResults = new Hashtable<String,String>();
		Hashtable<String,Float> mascotScores = new Hashtable<String,Float>();
		BufferedLineReader in = new BufferedLineReader(mascotFileName);
		String s;
		int prevScanNum = -1;
		while((s=in.readLine()) != null)
		{
			if(!s.startsWith("Elution"))
				continue;
			String[] token = s.split("\t");
			if(token.length != 5)
				continue;
			int scanNum = Integer.parseInt(token[0].substring(token[0].lastIndexOf(':')+2));
			if(scanNum == prevScanNum)
				continue;
			else
				prevScanNum = scanNum;
			int indexFileNum = token[0].lastIndexOf(".raw")-6;
			int fileNum = Integer.parseInt(token[0].substring(indexFileNum, indexFileNum+2));
			float score = Float.parseFloat(token[3]);
			mascotScores.put(fileNum+":"+scanNum, score);
			int c = Integer.parseInt(token[1]);
			if(c != charge)
				continue;
			if(score <= mascotThresholds.get(method+enzyme)[c-2])
				continue;
			String pep = token[2].substring(token[2].indexOf('.')+1, token[2].lastIndexOf('.'));
			mascotResults.put(fileNum+":"+scanNum, pep);
		}

		HashSet<String> mascotPepSet = new HashSet<String>();
		for(String pep : mascotResults.values())
		{
			StringBuffer seq = new StringBuffer();
			for(int i=0; i<pep.length(); i++)
			{
				char c = pep.charAt(i);
				if(c == 'Q')
					seq.append('K');
				else if(c == 'I')
					seq.append('L');
				else
					seq.append(c);
			}
			mascotPepSet.add(seq.toString());
		}
		
		Hashtable<String,String> msgfResults = new Hashtable<String,String>();
		Hashtable<String,Float> msgfScores = new Hashtable<String,Float>();
		
		File dir = new File(System.getProperty("user.home")+"/Research/Data/HeckWhole/Results0209");
		for(File f : dir.listFiles())
		{
			if(f.getName().endsWith("_"+db+".txt") && f.getName().contains(method) && f.getName().contains(enzyme))
			{
				prevScanNum = -1;
				in = new BufferedLineReader(f.getPath());
				while((s=in.readLine()) != null)
				{
					String[] token = s.split("\t");
					int scanNum = Integer.parseInt(token[1]);
					if(scanNum == prevScanNum)
						continue;
					else
						prevScanNum = scanNum;
					int fileNum = Integer.parseInt(f.getName().substring(f.getName().lastIndexOf("_"+db)-2, f.getName().lastIndexOf("_"+db)));
					int c = Integer.parseInt(token[4]);
					if(c != charge)
						continue;
					String pep = token[5].substring(token[5].indexOf('.')+1, token[5].lastIndexOf('.'));
					float specProb = Float.parseFloat(token[9]);
					float threshold = msgfThresholds.get(method+enzyme)[c-2];
					msgfScores.put(fileNum+":"+scanNum, specProb);
					if(specProb < threshold)
						msgfResults.put(fileNum+":"+scanNum, pep);
				}
			}
		}
		
		HashSet<String> msgfPepSet = new HashSet<String>();
		for(String pep : msgfResults.values())
		{
			StringBuffer seq = new StringBuffer();
			for(int i=0; i<pep.length(); i++)
			{
				char c = pep.charAt(i);
				if(c == 'Q')
					seq.append('K');
				else if(c == 'I')
					seq.append('L');
				else
					seq.append(c);
			}
			msgfPepSet.add(seq.toString());
		}		
//		System.out.print(method+"\t"+enzyme+"\t" + charge +"\t");
//		printVennDiagram(mascotResults, msgfResults);
		printVennDiagram(mascotPepSet, msgfPepSet);
		System.out.println();
//		String rawFileName;
//		if(enzyme.equalsIgnoreCase("Tryp"))
//			rawFileName = "090121_NM_Trypsin_";
//		else
//			rawFileName = "090128_NM_LysN_";
//		printSet1Exclusive(msgfResults, mascotResults, rawFileName, mascotScores, msgfScores);
	}

	private static void printSet1Exclusive(Hashtable<String,String> set1, Hashtable<String,String> set2, String rawFileName,
			Hashtable<String,Float> mascotScores, Hashtable<String,Float> msgfScores)
	{
		for(String s : set1.keySet())
		{
			if(set2.get(s) == null)	// set1 exclusive
			{
				String fileNum = s.split(":")[0];
				String scanNum = s.split(":")[1];
				Float score1 = mascotScores.get(s);
				Float score2 = msgfScores.get(s);
				System.out.println(rawFileName+fileNum+".RAW"+","+scanNum+","+set1.get(s)+","+score1+","+score2);
			}
		}
	}	
	
	private static void printVennDiagram(Hashtable<String,String> set1, Hashtable<String,String> set2)
	{
		int set1Only = 0;
		int set2Only = 0;
		int intersection = 0;
		int conflict = 0;
		
		for(String s : set1.keySet())
		{
			if(set2.get(s) != null)	// intersection
			{
				intersection++;
				String pep1 = set1.get(s);
				String pep2 = set2.get(s);
				boolean match = true;
				
				if(pep1.length() == pep2.length())
				{
					for(int i=0; i<pep1.length(); i++)
					{
						char c1 = pep1.charAt(i);
						char c2 = pep2.charAt(i);
						if(!(c1 == c2 || c1=='Q' && c2=='K' || c1=='I'&&c2=='L'
							|| c1=='K' && c2=='Q' || c1=='L'&&c2=='I'))
						{
							match = false;
							break;
						}
					}
				}
				else
					match = false;
				
				if(!match)
				{
//					System.out.println(s+" "+pep1+" "+pep2);
					conflict++;
				}
			}
			else
			{
				set1Only++;
			}
		}
		set2Only = set2.size()-intersection;
		System.out.print(set1Only+"\t"+set2Only+"\t"+intersection + "\t" + conflict);
	}
	
	private static void printVennDiagram(Set<String> set1, Set<String> set2)
	{
		int set1Only = 0;
		int set2Only = 0;
		int intersection = 0;
		
		for(String s : set1)
		{
			if(set2.contains(s))	// intersection
			{
				intersection++;
			}
			else
			{
				set1Only++;
//				System.out.println(s);
			}
		}
		set2Only = set2.size()-intersection;
		System.out.print(set1Only+"\t"+set2Only+"\t"+intersection);
	}
	
	public static void test() throws Exception
	{
		String method = "ETDTryp";
		String database = "_Target.txt";
		int charge = 2;
		
		File dir = new File(System.getProperty("user.home")+"/Research/Data/HeckWhole/Results0203");
		HashSet<String> pepset = new HashSet<String>();
		for(File f : dir.listFiles())
		{
			if(f.getName().endsWith(database) && f.getName().contains(method))
			{
				int prevScanNum = -1;
				BufferedLineReader in = new BufferedLineReader(f.getPath());
				String s;
				while((s=in.readLine()) != null)
				{
					String[] token = s.split("\t");
					int scanNum = Integer.parseInt(token[1]);
					if(scanNum == prevScanNum)
						continue;
					else
						prevScanNum = scanNum;
					int c = Integer.parseInt(token[4]);
					if(c != charge)
						continue;
					String pep = token[5].substring(token[5].indexOf('.')+1, token[5].lastIndexOf('.'));
					float specProb = Float.parseFloat(token[9]);
					float threshold = msgfThresholds.get(method)[c-2];
					if(specProb <= threshold)
						pepset.add(pep);
				}
			}
		}
		System.out.println(method+"\t"+charge+"\t"+pepset.size());
	}
	
	public static void vennDiagram() throws Exception
	{
		int charge = -1;
//		String method = "CID";
		String enzyme = "Tryp";
//		String enzyme = "LysN";
		
		Hashtable<String,String> cidIDTarget = new Hashtable<String,String>();
		Hashtable<String,String> cidIDDecoy = new Hashtable<String,String>();
		File dir = new File(System.getProperty("user.home")+"/Research/Data/HeckRevision/MSGFDB0720_AAFreq");
		HashSet<String> cid = new HashSet<String>();
		HashSet<String> etd = new HashSet<String>();
		for(File f : dir.listFiles())
		{
			if(f.getName().endsWith("_0.txt") && f.getName().contains(enzyme) && 
					(f.getName().contains("CID")))
			{
				BufferedLineReader in = new BufferedLineReader(f.getPath());
				String s;
				int prevScanNum = -1;
				while((s=in.readLine()) != null)
				{
					if(s.startsWith("#"))
						continue;
					String[] token = s.split("\t");
					int scanNum = Integer.parseInt(token[1]);
					if(scanNum == prevScanNum)
						continue;
					else
						prevScanNum = scanNum;
					int c = Integer.parseInt(token[4]);
					if(c > 4)
						c = 4;
					if(charge > 0 && c != charge)
						continue;
					String pep = token[5].substring(token[5].indexOf('.')+1, token[5].lastIndexOf('.'))+":"+c;
					float specProb = Float.parseFloat(token[10]);
					if(specProb < msgfThresholds.get("CID"+enzyme)[c-2])
					{
						cid.add(pep);
						cidIDTarget.put(token[0]+":"+token[1], pep);
					}
				}
			}
			else if(f.getName().endsWith("_1.txt") && f.getName().contains(enzyme)
					&& f.getName().contains("CID"))
			{
				BufferedLineReader in = new BufferedLineReader(f.getPath());
				String s;
				int prevScanNum = -1;
				while((s=in.readLine()) != null)
				{
					if(s.startsWith("#"))
						continue;
					String[] token = s.split("\t");
					int scanNum = Integer.parseInt(token[1]);
					if(scanNum == prevScanNum)
						continue;
					else
						prevScanNum = scanNum;
					int c = Integer.parseInt(token[4]);
					if(c > 4)
						c = 4;
					if(charge > 0 && c != charge)
						continue;
					
					String pep = token[5].substring(token[5].indexOf('.')+1, token[5].lastIndexOf('.'))+":"+c;
					float specProb = Float.parseFloat(token[10]);
					if(specProb < msgfThresholds.get("CID"+enzyme)[c-2])
					{
						cidIDDecoy.put(token[0]+":"+token[1], pep);
					}
				}
			}
		}
		
		int numCIDOnlyTarget = 0;
		HashSet<String> cidOnlyTarget = new HashSet<String>();
		int numCIDOnlyDecoy = 0;
		HashSet<String> cidOnlyDecoy = new HashSet<String>();
		int numETDOnlyTarget = 0;
		HashSet<String> etdOnlyTarget = new HashSet<String>();
		int numETDOnlyDecoy = 0;
		HashSet<String> etdOnlyDecoy = new HashSet<String>();
		int numBothTarget = 0;
		HashSet<String> bothTarget = new HashSet<String>();
		int numBothDecoy = 0;
		HashSet<String> conflictTarget = new HashSet<String>();
		HashSet<String> bothDecoy = new HashSet<String>();
		HashSet<String> conflictDecoy = new HashSet<String>();
		int numConflictTarget = 0;
		int numConflictDecoy = 0;
		
		for(File f : dir.listFiles())
		{
			if(f.getName().endsWith("_0.txt") && f.getName().contains(enzyme) && 
					(f.getName().contains("ETD")))
			{
				BufferedLineReader in = new BufferedLineReader(f.getPath());
				String s;
				int prevScanNum = -1;
				while((s=in.readLine()) != null)
				{
					if(s.startsWith("#"))
						continue;
					String[] token = s.split("\t");
					int scanNum = Integer.parseInt(token[1]);
					if(scanNum == prevScanNum)
						continue;
					else
						prevScanNum = scanNum;
					int c = Integer.parseInt(token[4]);
					if(c > 4)
						c = 4;
					if(charge > 0 && c != charge)
						continue;
					
					String pep = token[5].substring(token[5].indexOf('.')+1, token[5].lastIndexOf('.'))+":"+c;
					float specProb = Float.parseFloat(token[10]);
					if(specProb < msgfThresholds.get("ETD"+enzyme)[c-2])
					{
						etd.add(pep);
						String cidPep = cidIDTarget.get(token[0]+":"+(scanNum-1));
						if(cidPep == null)
						{
							numETDOnlyTarget++;
							etdOnlyTarget.add(pep);
						}
						else
						{
							if(cidPep.equalsIgnoreCase(pep))	// match
							{
								numBothTarget++;
								bothTarget.add(pep);
							}
							else							// conflict
							{
								numConflictTarget++;
								conflictTarget.add(cidPep+":"+pep);
							}
							cidIDTarget.remove(token[0]+":"+(scanNum-1));
						}
					}
				}
			}
			else if(f.getName().endsWith("_1.txt") && f.getName().contains(enzyme)
					&& f.getName().contains("ETD"))
			{
				BufferedLineReader in = new BufferedLineReader(f.getPath());
				String s;
				int prevScanNum = -1;
				while((s=in.readLine()) != null)
				{
					if(s.startsWith("#"))
						continue;
					String[] token = s.split("\t");
					int scanNum = Integer.parseInt(token[1]);
					if(scanNum == prevScanNum)
						continue;
					else
						prevScanNum = scanNum;
					int c = Integer.parseInt(token[4]);
					if(c > 4)
						c = 4;
					if(charge > 0 && c != charge)
						continue;
					
					String pep = token[5].substring(token[5].indexOf('.')+1, token[5].lastIndexOf('.'))+":"+c;
					float specProb = Float.parseFloat(token[10]);
					if(specProb < msgfThresholds.get("ETD"+enzyme)[c-2])
					{
						String cidPep = cidIDDecoy.get(token[0]+":"+(scanNum-1));
						if(cidPep == null)
						{
							numETDOnlyDecoy++;
							etdOnlyDecoy.add(pep);
						}
						else
						{
							if(cidPep.equalsIgnoreCase(pep))	// match
							{
								numBothDecoy++;
								bothDecoy.add(pep);
							}
							else							// conflict
							{
								numConflictDecoy++;
								conflictDecoy.add(cidPep+":"+pep);
							}
							cidIDDecoy.remove(token[0]+":"+(scanNum-1));
						}
					}
				}				
			}
		}
		
		numCIDOnlyTarget = cidIDTarget.size();
		for(String pep : cidIDTarget.values())
			cidOnlyTarget.add(pep);
		numCIDOnlyDecoy = cidIDDecoy.size();
		for(String pep : cidIDDecoy.values())
			cidOnlyDecoy.add(pep);
				
		int numPepCIDOnlyTarget = cidOnlyTarget.size();
		for(String pep : cidOnlyTarget)
			if(bothTarget.contains(pep))
				numPepCIDOnlyTarget--;

		int numPepETDOnlyTarget = etdOnlyTarget.size();
		for(String pep : etdOnlyTarget)
			if(bothTarget.contains(pep))
				numPepETDOnlyTarget--;

		int numPepCIDOnlyDecoy = cidOnlyDecoy.size();
		for(String pep : cidOnlyDecoy)
			if(bothDecoy.contains(pep))
				numPepCIDOnlyDecoy--;

		int numPepETDOnlyDecoy = etdOnlyDecoy.size();
		for(String pep : etdOnlyDecoy)
			if(bothDecoy.contains(pep))
				numPepETDOnlyDecoy--;
		
		System.out.println(enzyme);
		System.out.println("*****Target");
		System.out.println("NumCIDOnly: " + numCIDOnlyTarget + "\t" + numPepCIDOnlyTarget);
		System.out.println("NumETDOnly: " + numETDOnlyTarget+ "\t" + numPepETDOnlyTarget);
		System.out.println("NumIntersection: " + (numBothTarget+numConflictTarget)+ "\t" + bothTarget.size());
		System.out.println("NumConflict: " + numConflictTarget+"\t"+conflictTarget.size());
		System.out.println("*****Decoy");
		System.out.println("NumCIDOnly: " + numCIDOnlyDecoy+ "\t" + numPepCIDOnlyDecoy);
		System.out.println("NumETDOnly: " + numETDOnlyDecoy+ "\t" + numPepETDOnlyDecoy);
		System.out.println("NumIntersection: " + (numBothDecoy+numConflictDecoy)+ "\t" + bothDecoy.size());
		System.out.println("NumConflict: " + numConflictDecoy + "\t"+conflictDecoy.size());
		System.out.println(cid.size()+"\t"+etd.size());
	}
	
	public static void extractGoodSpectraFromSum() throws Exception
	{
		String[] enzyme = {"Tryp", "LysN"};
		for(String e : enzyme)
			extractGoodSpectraFromSum(e, 
					System.getProperty("user.home")+"/Research/Data/HeckWhole/AnnotatedSpectra/annotatedHeckWholeSum_CID_"+e+".mgf", 
					System.getProperty("user.home")+"/Research/Data/HeckWhole/AnnotatedSpectra/annotatedHeckWholeSum_ETD_"+e+".mgf");
	}
	
	public static void rankOnePeakTest() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/HeckWhole/AnnotatedSpectra/annotatedHeckWhole_ETD_Tryp.mgf";
		int charge = 2;
		Tolerance tolerance = new Tolerance(0.5f);
		
		IonType[] ionType = {IonType.A, IonType.B, IonType.C, IonType.Y, IonType.Z,
				IonType.getIonType("b2"), IonType.getIonType("c2"),
				IonType.getIonType("y2"), IonType.getIonType("z2"),
				IonType.getIonType("c+H"), IonType.getIonType("z+H"),
				IonType.getIonType("c-H"), IonType.getIonType("z-H"),
				IonType.getIonType("b+H"), IonType.getIonType("y+H"),
				IonType.getIonType("b-H"), IonType.getIonType("y-H")
				};
		
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			if(spec.getCharge() != charge)
				continue;
			
			System.out.println("*****"+spec.getAnnotationStr()+"\t"+spec.getCharge());
			spec.setActivationMethod(ActivationMethod.ETD);
			spec.setRanksOfPeaks();
			Peak rankOnePeak = null;
			for(Peak p : spec)
				if(p.getRank() == 1)
					rankOnePeak = p;
			
			
			boolean isExplained = false;
			spec.setRanksOfPeaks();
			Spectrum newSpec = spec.getCloneWithoutPeakList();
			newSpec.add(rankOnePeak);
			SpectrumAnnotator annotator = new SpectrumAnnotator(newSpec, spec.getAnnotation());
			for(IonType ion : ionType)
			{
				if(annotator.getPeakListOfIon(ion, tolerance).size() > 0)
				{
					System.out.println(ion.getName());
					isExplained = true;
					break;
				}
			}
			
			
			if(!isExplained)
			{
				for(int c=1; c<=spec.getCharge(); c++)
				{
					float mass = (spec.getParentMass()+c*(float)(Composition.NEUTRON))/c;
					System.out.println("NE " + c + " " + (rankOnePeak.getMz()-mass) + " " + rankOnePeak.getMz() + " " + mass);
				}		
			}
			
		}
	}
	
	public static void precursorTest() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/HeckWhole/AnnotatedSpectra/annotatedHeckWhole_ETD_LysN.mgf";
		int charge = 3;
		Tolerance tolerance = new Tolerance(0.3f);
		
		float[] offsets = {
				0,
				(float)Composition.ISOTOPE,
				(float)Composition.ISOTOPE2,
				-(float)Composition.H2O,
				-(float)Composition.H2O*2,
				-(float)Composition.NH3,
				-(float)Composition.NH3*2
		};
		
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			spec.setRanksOfPeaks();
			if(spec.getCharge() != charge)
				continue;
			System.out.println(spec.getAnnotationStr()+"\t"+spec.getCharge());
			float precursorMz = spec.getPrecursorPeak().getMz();
			// remove charge reduced precursor masses
			float neutralParentMass = (precursorMz-(float)Composition.NEUTRON)*spec.getCharge();
			int startCharge = 2;
			for(int c=startCharge; c<=spec.getCharge(); c++)
			{
				for(float off : offsets)
				{
					float mass = (neutralParentMass+c*(float)Composition.NEUTRON+off)/c;
					for(Peak p : spec.getPeakListByMass(mass, tolerance))
					{
						if(p.getRank() <= 10)
						{
							System.out.println(c+"\t"+off+"\t"+p.getRank());
							p.setIntensity(0);
							spec.setRanksOfPeaks();
						}
					}
				}
			}			
		}
	}
	public static void countNumSpectra() throws Exception
	{
		File dir = new File("/home/sangtaekim/Research/Data/HeckWhole/Spectra");
		int[][] numSpec = new int[2][100];	// numSpec[0]: Trypsin, numSpec[1]: LysN
		for(File f : dir.listFiles())
		{
			if(f.getName().endsWith(".mzXML") || f.getName().endsWith(".mgf"))
			{
				int frag = 0;
				if(f.getName().contains("LysN"))
					frag = 1;
				Iterator<Spectrum> itr = null;
				if(f.getName().endsWith(".mzXML"))
					itr = new MzXMLSpectraIterator(f.getPath());
				else
					itr = new edu.ucsd.msjava.msutil.SpectraIterator(f.getPath(), new MgfSpectrumParser());
				int prevScanNum = -1;
				float prevPrecursorMz = 0;
				while(itr.hasNext())
				{
					Spectrum spec = itr.next();
					if(spec.getScanNum() == prevScanNum+1 && Math.abs(prevPrecursorMz-spec.getPrecursorPeak().getMz()) < 0.01f)	// paired
					{
						numSpec[frag][spec.getCharge()]++;
					}
					else
					{
						prevScanNum = spec.getScanNum();
						prevPrecursorMz = spec.getPrecursorPeak().getMz();
					}
				}
			}
		}
		for(int i=0; i<=1; i++)
		{
			String method;
			if(i == 0)
				method = "Trypsin";
			else
				method = "LysN";
			System.out.println(method);
			for(int j=1; j<100; j++)
			{
				if(numSpec[i][j] > 0)
					System.out.println(j+"\t"+numSpec[i][j]);
			}
		}
		System.out.println("Done");
	}
	
	public static void convertMzXMLIntoMgf() throws Exception
	{
		String dirName = System.getProperty("user.home")+"/Research/Data/HeckWhole/Spectra";
		File dir = new File(dirName);
		
		ActivationMethod[] methods = {ActivationMethod.CID, ActivationMethod.ETD};
		String[] enzymes = {"Tryp", "LysN"};
		
		for(String enzyme : enzymes)
		{
			for(ActivationMethod method : methods)
			{
				for(File f : dir.listFiles())
				{
					if(f.getName().contains(enzyme) && f.getName().endsWith(".mzXML"))
					{
						System.out.print(f.getName()+" " + method + " " + enzyme + "...");
						String mgfFileName = f.getParent()+File.separator+f.getName().substring(0, f.getName().lastIndexOf('.'))+"_"+method+".mgf";
						MzXMLToMgfConverter.convert(f, new File(mgfFileName), 0, 10000, method, 2, 2);
					}
				}
			}
		}
	}

	public static void extractGoodSpectra() throws Exception
	{
		String[] methods = {"CID", "ETD"};
		String[] enzymes = {"Tryp","LysN"};
		for(String method : methods)
		{
			for(String enzyme : enzymes)
			{
				extractGoodSpectra(method, enzyme, System.getProperty("user.home")+"/Research/Data/HeckRevision/AnnotatedSpectra/"+method+"_"+enzyme+"_Confident.mgf");
			}
		}
		System.out.println("Done");
	}
	
	public static void countEnzymeCleavage(String method, String enzyme) throws Exception
	{
		File dir = new File(System.getProperty("user.home")+"/Research/Data/HeckRevision");
		File resDir = new File(dir.getPath()+"/MSGFDB0720_AAFreq");
		
		String metID = method+enzyme;
		int numID = 0;
		int numPrecedingKR = 0;
		int numEndingKR = 0;
		for(File f : resDir.listFiles())
		{
			if(f.getName().endsWith("_0.txt") && f.getName().contains(metID))
			{
				BufferedLineReader in = new BufferedLineReader(f.getPath());
				String s;
				int prevScanNum = -1;
				while((s=in.readLine()) != null)
				{
					if(s.startsWith("#"))
						continue;
					String[] token = s.split("\t");
					if(token.length != 11)
						assert(false);
					
					if(!token[2].equalsIgnoreCase(method))
						continue;
					int scanNum = Integer.parseInt(token[1]);
					if(scanNum == prevScanNum)
						continue;
					else
						prevScanNum = scanNum;

					int charge = Integer.parseInt(token[4]);
					int chargeIndex = charge-2;
					if(chargeIndex > 2)
						chargeIndex = 2;
					
					float specProb = Float.parseFloat(token[10]);
					if(specProb >= msgfThresholds.get(metID)[chargeIndex])
						continue;
					String annotation = token[5];
					numID++;
					if(annotation.charAt(0) == 'K' || annotation.charAt(0) == 'R')
						++numPrecedingKR;
					char lastAA = annotation.charAt(annotation.lastIndexOf('.')-1);
					if(lastAA == 'K' || lastAA == 'R')
						++numEndingKR;
				}
			}
		}
		System.out.println(numID+"\t"+numPrecedingKR/(float)numID+"\t"+numEndingKR/(float)numID);
	}
	
	public static void extractGoodSpectra(String method, String enzyme, String outFileName) throws Exception
	{
		System.out.println(method+enzyme+" "+outFileName);
		File dir = new File(System.getProperty("user.home")+"/Research/Data/HeckRevision");
		File resDir = new File(dir.getPath()+"/MSGFDB0720_AAFreq");
		File specDir = new File(System.getProperty("user.home")+"/Research/Data/HeckWhole/Spectra");
		
		String metID = method+enzyme;
		
		Hashtable<String, SpectrumAccessorBySpecIndex> specTable = new Hashtable<String, SpectrumAccessorBySpecIndex>(); 
		for(File f : specDir.listFiles())
		{
			if(f.getName().contains(enzyme))
			{
				if(f.getName().endsWith("mzXML"))
					specTable.put(f.getName(), new MzXMLSpectraMap(f.getPath()));
			}
		}
		
		Hashtable<String, Float> pepScoreTable = new Hashtable<String, Float>();	// stores the best spec prob of peptides
		Hashtable<String, String> pepSpecTable = new Hashtable<String, String>();	// stores the best spec of peptides
		
		for(File f : resDir.listFiles())
		{
			if(f.getName().endsWith("_0.txt") && f.getName().contains(metID))
			{
				BufferedLineReader in = new BufferedLineReader(f.getPath());
				String s;
				int prevScanNum = -1;
				while((s=in.readLine()) != null)
				{
					if(s.startsWith("#"))
						continue;
					String[] token = s.split("\t");
					if(token.length != 11)
						assert(false);
					
					if(!token[2].equalsIgnoreCase(method))
						continue;
					int scanNum = Integer.parseInt(token[1]);
					if(scanNum == prevScanNum)
						continue;
					else
						prevScanNum = scanNum;

					int charge = Integer.parseInt(token[4]);
					int chargeIndex = charge-2;
					if(chargeIndex > 2)
						chargeIndex = 2;
//					if(charge != 2 && charge != 3)
//						continue;
					
					float specProb = Float.parseFloat(token[10]);
					if(specProb >= msgfThresholds.get(metID)[chargeIndex])
						continue;
					String peptide = token[5].substring(token[5].indexOf('.')+1, token[5].lastIndexOf('.'));
//					String peptide = token[5];
					Float prevBestScore = pepScoreTable.get(peptide+"|"+chargeIndex);	// peptide:charge
					if(prevBestScore == null || specProb < prevBestScore)
					{
						pepScoreTable.put(peptide+"|"+chargeIndex, specProb);
						pepSpecTable.put(peptide+"|"+chargeIndex, token[0]+":"+token[1]+":"+token[7]+":"+token[8]);
//						pepSpecTable.put(peptide, token[0]+":"+(Integer.parseInt(token[1])-1));
					}
				}
			}
		}
		
		PrintStream mgfOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(outFileName)));
		ArrayList<String> pepSet = new ArrayList<String>(pepSpecTable.keySet());
		for(String pep : pepSet)
		{
			String specInfo = pepSpecTable.get(pep);
			String[] token = specInfo.split(":");
			assert(token.length == 4): specInfo;
			String fileName = token[0];
			int scanNum = Integer.parseInt(token[1]);
			int msgfScore = Integer.parseInt(token[2]);
			int peptideScore = (int)Float.parseFloat(token[3]);
			SpectrumAccessorBySpecIndex map = specTable.get(fileName);
			assert(map != null);
			Spectrum spec = map.getSpectrumBySpecIndex(scanNum);
			assert(spec != null);
			Float specProb = pepScoreTable.get(pep);
			assert(specProb != null);
			spec.setTitle(fileName+":"+scanNum+":"+specProb+":"+msgfScore+":"+peptideScore);
//			String pepSeq = pep.substring(pep.indexOf('.')+1, pep.lastIndexOf('.'));
			String pepSeq = pep.substring(0, pep.lastIndexOf('|'));
			spec.setAnnotation(new Peptide(pepSeq, AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys()));
			spec.outputMgf(mgfOut);
		}
		mgfOut.flush();
		mgfOut.close();
	}	
	public static void extractGoodSpectraFromMascot(String method, String enzyme, String outFileName) throws Exception
	{
		File dir = new File(System.getProperty("user.home")+"/Research/Data/HeckWhole/");
		File resDir = new File(dir.getPath()+"/ResultsShabaz");
		File specDir = new File(dir.getPath()+"/Spectra");
		
		String metID = method+enzyme;
		
		Hashtable<String, SpectrumAccessorBySpecIndex> specTable = new Hashtable<String, SpectrumAccessorBySpecIndex>(); 
		for(File f : specDir.listFiles())
		{
			if(f.getName().contains(enzyme))
			{
				if(f.getName().endsWith("mzXML"))
					specTable.put(f.getName(), new MzXMLSpectraMap(f.getPath()));
				else if(f.getName().endsWith("mgf"))
					specTable.put(f.getName(), new SpectraMap(f.getPath(), new MgfSpectrumParser()));
			}
		}
		
		Hashtable<String, Float> pepScoreTable = new Hashtable<String, Float>();	// stores the best spec prob of peptides
		Hashtable<String, String> pepSpecTable = new Hashtable<String, String>();	// stores the best spec of peptides
		
		for(File f : resDir.listFiles())
		{
			if(f.getName().endsWith("_Target.txt") && f.getName().contains(method) && f.getName().contains(enzyme))
			{
				BufferedLineReader in = new BufferedLineReader(f.getPath());
				String s;
				String prevScanNum = "";
				while((s=in.readLine()) != null)
				{
					if(s.startsWith("#"))
						continue;
					String[] token = s.split("\t");
					if(token.length != 5)
						assert(false);
					
					String scanNum = token[0];
					if(scanNum.equalsIgnoreCase(prevScanNum))
						continue;
					else
						prevScanNum = scanNum;

					int charge = Integer.parseInt(token[1]);
					if(charge > 4)
						charge = 4;
					float mascotScore = Float.parseFloat(token[9]);
					if(mascotScore >= mascotThresholds.get(metID)[charge-2])
						continue;
					String peptide = token[2].substring(token[2].indexOf('.')+1, token[2].lastIndexOf('.'));
//					String peptide = token[5];
					Float prevBestScore = pepScoreTable.get(peptide+"|"+token[1]);	// peptide:charge
					if(prevBestScore == null || mascotScore < prevBestScore)
					{
						pepScoreTable.put(peptide+"|"+token[1], mascotScore);
						pepSpecTable.put(peptide+"|"+token[1], token[0]);
//						pepSpecTable.put(peptide, token[0]+":"+(Integer.parseInt(token[1])-1));
					}
				}
			}
		}
		
		PrintStream mgfOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(outFileName)));
		ArrayList<String> pepSet = new ArrayList<String>(pepSpecTable.keySet());
		for(String pep : pepSet)
		{
			String specInfo = pepSpecTable.get(pep);
			String fileName = specInfo.substring(specInfo.indexOf("_nm_")+4, specInfo.lastIndexOf("raw")-5);
			int scanNum = Integer.parseInt(specInfo.substring(specInfo.lastIndexOf("ScanNumber:")+12));
			SpectrumAccessorBySpecIndex map = specTable.get(fileName);
			assert(map != null);
			Spectrum spec = map.getSpectrumBySpecIndex(scanNum);
			assert(spec != null);
			Float mascotScore = pepScoreTable.get(pep);
			assert(mascotScore != null);
			spec.setTitle(specInfo+":"+mascotScore);
//			String pepSeq = pep.substring(pep.indexOf('.')+1, pep.lastIndexOf('.'));
			String pepSeq = pep.substring(0, pep.lastIndexOf('|'));
			spec.setAnnotation(new Peptide(pepSeq, AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys()));
			spec.outputMgf(mgfOut);
		}
		mgfOut.flush();
		mgfOut.close();
	}	
	
	public static void extractGoodSpectraFromSum(String enzyme, String cidOutFileName, String etdOutFileName) throws Exception
	{
		File dir = new File(System.getProperty("user.home")+"/Research/Data/HeckWhole");
		File resDir = new File(dir.getPath()+"/Results0209");
		File specDir = new File(dir.getPath()+"/Spectra");
		
		Hashtable<String, SpectrumAccessorBySpecIndex> specTable = new Hashtable<String, SpectrumAccessorBySpecIndex>(); 
		for(File f : specDir.listFiles())
		{
			if(f.getName().endsWith("mzXML"))
				specTable.put(f.getName(), new MzXMLSpectraMap(f.getPath()));
			else if(f.getName().endsWith("mgf"))
				specTable.put(f.getName(), new SpectraMap(f.getPath(), new MgfSpectrumParser()));
		}
		
		Hashtable<String, Float> pepScoreTable = new Hashtable<String, Float>();	// stores the best spec prob of peptides
		Hashtable<String, String> pepSpecTable = new Hashtable<String, String>();	// stores the best spec of peptides
		
		for(File f : resDir.listFiles())
		{
			if(f.getName().contains("Sum") && f.getName().contains(enzyme) && f.getName().endsWith("_Target.txt"))
			{
				BufferedLineReader in = new BufferedLineReader(f.getPath());
				String s;
				String prevScanNum = "";
				while((s=in.readLine()) != null)
				{
					if(s.startsWith("#"))
						continue;
					String[] token = s.split("\t");
					if(token.length != 10)
						assert(false);
					
					if(!token[2].equalsIgnoreCase("CID/ETD"))
						continue;

					int charge = Integer.parseInt(token[4]);
					if(charge != 2 && charge != 3)
						continue;
					
					String scanNum = token[1];
					if(scanNum.equalsIgnoreCase(prevScanNum))
						continue;
					else
						prevScanNum = scanNum;
					
					float specProb = Float.parseFloat(token[9]);
					if(specProb >= msgfThresholds.get("Sum"+enzyme)[charge-2])
						continue;
					
					String peptide = token[5];
					peptide = peptide.substring(peptide.indexOf('.')+1, peptide.lastIndexOf('.'));
					Float prevBestScore = pepScoreTable.get(peptide+"|"+token[4]);	// peptide:charge
					if(prevBestScore == null || specProb < prevBestScore)
					{
						pepScoreTable.put(peptide+"|"+token[4], specProb);
						pepSpecTable.put(peptide+"|"+token[4], token[0]+":"+token[1]+":"+token[7]+":"+token[8]);
					}
				}
			}
		}
		
		PrintStream cidOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(cidOutFileName)));
		PrintStream etdOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(etdOutFileName)));
//		SpectraContainer containerCID = new SpectraContainer();
//		SpectraContainer containerETD = new SpectraContainer();
		ArrayList<String> pepSet = new ArrayList<String>(pepSpecTable.keySet());
		int num = 0;
		for(String pep : pepSet)
		{
//			System.out.println(num++);
			String specInfo = pepSpecTable.get(pep);
			String[] token = specInfo.split(":");
			assert(token.length == 4): specInfo;
			String fileName = token[0];
			int scanNumCID = Integer.parseInt(token[1].substring(0, token[1].indexOf('-')));
			int scanNumETD = scanNumCID+1;
			int msgfScore = Integer.parseInt(token[2]);
			int peptideScore = (int)Float.parseFloat(token[3]);
			
			SpectrumAccessorBySpecIndex map = specTable.get(fileName);
			assert(map != null);
			
			Spectrum specCID = map.getSpectrumBySpecIndex(scanNumCID);
			Spectrum specETD = map.getSpectrumBySpecIndex(scanNumETD);
			assert(specCID != null && specETD != null);
			Float specProb = pepScoreTable.get(pep);
			assert(specProb != null);
			specCID.setTitle(fileName+":CID:"+scanNumCID+":"+specProb+":"+msgfScore+":"+peptideScore);
			specETD.setTitle(fileName+":ETD:"+scanNumETD+":"+specProb+":"+msgfScore+":"+peptideScore);
			String pepSeq = pep.substring(0, pep.lastIndexOf('|'));
			Peptide annotation = new Peptide(pepSeq, AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys()); 
			specCID.setAnnotation(annotation);
			specETD.setAnnotation(annotation);
			specCID.outputMgf(cidOut);
			specETD.outputMgf(etdOut);
		}
		
		cidOut.flush();
		cidOut.close();
		etdOut.flush();
		etdOut.close();
		System.out.println("Done: " + enzyme + " " + cidOutFileName + " " + etdOutFileName);
	}	
}
