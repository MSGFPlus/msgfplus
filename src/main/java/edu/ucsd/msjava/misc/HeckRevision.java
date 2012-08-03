package edu.ucsd.msjava.misc;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;

import edu.ucsd.msjava.msgf.AminoAcidGraph;
import edu.ucsd.msjava.msgf.GeneratingFunction;
import edu.ucsd.msjava.msgf.Histogram;
import edu.ucsd.msjava.msgf.ScoredSpectrum;
import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.NewScorerFactory;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.parser.MascotParser;
import edu.ucsd.msjava.parser.MgfSpectrumParser;
import edu.ucsd.msjava.parser.MzXMLSpectraIterator;
import edu.ucsd.msjava.parser.MzXMLSpectraMap;

public class HeckRevision {
	public static void main(String argv[]) throws Exception
	{
//		analyzeMascotResults();
//		splitTargetDecoy();
//		splitMascotResult();
//		processMSGDResult();
//		compareMzXMLAndMgf();
//		analyzeMascotIDs();
//		analyzeMSGDIDs();
//		mergeSpectralPairs();
//		makeOMSSAScript();
//		makeOMSSAScriptMerged();
//		makeSampleSpectra();
//		processPercolatorResults();
//		makeRescoringScript("ETD", "Tryp");
//		specProbTest();
//		drawRankDist();
	}
	
	public static void drawRankDist() throws Exception 
	{
		HashMap<String,String> ionMap = new HashMap<String,String>();
		ionMap.put("S,1,19", "y");
		ionMap.put("P,1,1", "b");
		ionMap.put("S,1,20", "y+n");
		ionMap.put("P,1,-17", "b-H2O");
		ionMap.put("P,1,-16", "b-NH3");
		ionMap.put("P,1,2", "b+n");
		ionMap.put("S,1,1", "y-H2O");
		ionMap.put("S,1,2", "y-NH3");
		ionMap.put("S,2,10", "y2");
		ionMap.put("P,2,1", "b2");
		ionMap.put("S,2,11", "y2+n");
		ionMap.put("S,3,7", "y3");
		ionMap.put("S,1,5", "z+H2");
		ionMap.put("S,1,4", "z+n(z+H)");
		ionMap.put("S,1,3", "z");
		ionMap.put("P,1,17", "c-H");
		ionMap.put("P,1,18", "c");
		ionMap.put("P,1,19", "c+n(c+H)");
		ionMap.put("P,1,20", "c+H2");
		ionMap.put("S,2,2", "z2");
		ionMap.put("P,2,10", "c2");
		ionMap.put("S,2,3", "z2+H");
//		ionMap.put("", "");
//		ionMap.put("", "");
//		ionMap.put("", "");
//		ionMap.put("", "");
		
		String activation = "ETD";//ETD
		String enzyme = "Tryp";//LysN
		for(int charge = 2;charge<5;charge++){
			String filename = "/home/sangtaekim/Research/Data/HeckRevision/AnnotatedSpectra/plots/" + activation + enzyme +
			"Charge"+charge+".txt";
			String title = activation+"-"+enzyme + " charge " + charge;
			String outfilename = filename + ".m";
			BufferedLineReader in = new BufferedLineReader(filename);
			BufferedWriter out = new BufferedWriter(new FileWriter(outfilename));
			String s;
			boolean read = false;
			int numPartition = Integer.MAX_VALUE;
			int n = 1;
			ArrayList<String> ionNames = new ArrayList<String>();
			out.write("close;\nclear;\n");
			while((s=in.readLine())!=null){

				if(s.startsWith("#RankDistributions")){
					numPartition = Integer.parseInt(s.split("\t")[1]);
					read = true;
					n=1;
				}

				if(s.startsWith("#ErrorDistributions")) read = false;

				if(read){
					if(numPartition == 0){
						String[] token = s.split("\t");
						if(token[0].startsWith("noise")) continue;
						ionNames.add(token[0]);
						out.write("ionDist(:," +n + ")=[");
						System.out.print("ionDist(:," +(n++) + ")=[");
						for(int i=1; i< token.length-2;i++){
							out.write(token[i] + " ");
							System.out.print(token[i] + " ");
						}

						out.write("];\n");
						System.out.println("];");
					}


					if(s.startsWith("Partition")) numPartition--;

				}

			}
			ionNames.add("Unexplained");
			out.write("\nionDist(:," + n + ")= 1-sum(ionDist');\n");
			out.write("figure1 = figure('XVisual',...\n'0x24 (TrueColor, depth24, RGB mask 0xff0000 0xff00 0x00ff)');\n");
			out.write("axes1 = axes('Parent',figure1,'FontSize',20);\n");
			out.write("box('on');hold('all');\n");
			out.write("plot1 = plot(ionDist);\n");

			for(int i=0;i<ionNames.size();i++){
				String name = ionNames.get(i).replace('_',',');
				assert(name.equalsIgnoreCase("Unexplained") || ionMap.get(name) != null): name;
				if(i<ionNames.size()/2)
					out.write("set(plot1("+(i+1)+"),'DisplayName','"+ionNames.get(i).replace('_',',')+"','LineWidth',1);\n");
				else if(i<ionNames.size()-1)
					out.write("set(plot1("+(i+1)+"),'DisplayName','"+ionNames.get(i).replace('_',',')+"','Marker','*','LineStyle','none');\n");
				else
					out.write("set(plot1("+(i+1)+"),'DisplayName','"+ionNames.get(i).replace('_',',')+"','LineWidth',1,'Color',[0 0 0]);\n");
			}
			out.write("xlabel('Rank','FontSize',20);\n");
			out.write("ylabel('Probability','FontSize',20);\n");
			out.write("title({'"+title+"'},'FontSize',20,'FontName','helvetica');\n");
			out.write("legend(axes1,'show');\n");
			out.write("xlim([1 100]);\nylim([0 1])\n");
			out.close();
			in.close();
		}
	}

	public static void makeRescoringScript(String method, String enzyme) throws Exception
	{
		File inputDir = new File("/home/sangtaekim/Research/Data/HeckRevision/MSGFDB0720");
		File outputDir = new File("/home/sangtaekim/Research/Data/HeckRevision/MSGFDB0720_AAFreq");
		
		System.out.println("#!/bin/bash");
		for(File f : inputDir.listFiles())
		{
			String fileName = f.getName();
			if(fileName.startsWith("MSGFDB") && fileName.contains(method) && fileName.contains(enzyme))
			{
				System.out.println("java -jar /home/sangtaekim/Research/ToolDistribution/MSGF.jar " +
						"-i " + f.getPath() + " " +
						"-d /home/sangtaekim/Research/Data/HeckWhole/Spectra/ " +
						" -o " + outputDir.getPath() + File.separator + f.getName() + " " +  
						"-m " + (method.equalsIgnoreCase("ETD") ? 1 : 0) + " " +
						"-e " + (enzyme.equalsIgnoreCase("LysN") ? 4 : 1) + " " +
						"-aaSet /home/sangtaekim/Research/Data/HeckRevision/MSGFDB0720_AAFreq/AAMasses.txt");
			}
		}
	}
	
	public static void processPercolatorResults() throws Exception
	{
		File dirPercolator = new File("/home/sangtaekim/Research/Data/HeckRevision/Percolator");
		File dirMascot = new File("/home/sangtaekim/Research/Data/HeckRevision/Mascot23");
		
		for(File f : dirPercolator.listFiles())
		{
			String fileName = f.getName();
			if(!fileName.endsWith(".pop"))
				continue;
			System.out.println(fileName);
			String name = fileName.substring(0, fileName.indexOf('.'));
			String mascotResultFileName = dirMascot.getPath() + File.separator + name + ".dat";
			HashMap<Integer,Integer> map = MascotParser.getQueryNumChargeMap(mascotResultFileName);
			
			String db;
			if(fileName.contains("target"))
				db = "target";
			else
				db = "decoy";
			String outputFileName = dirPercolator + File.separator + name + "." + db + ".pout";
			PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFileName)));

			BufferedLineReader in = new BufferedLineReader(f.getPath());
			out.println("Charge\t"+in.readLine());	// header
			String s;
			while((s=in.readLine()) != null)
			{
				String[] token = s.split("\t");
				assert(token.length >= 6): s;
				int queryNum = Integer.parseInt(token[0].substring(token[0].indexOf(':')+1, token[0].indexOf(';')));
				int charge = map.get(queryNum);
				out.println(charge+"\t"+s);
			}
			in.close();
			out.close();
		}
		System.out.println("Done");
	}
	
	public static void makeSampleSpectra() throws Exception
	{
		String specFile = "/home/sangtaekim/Research/Data/HeckWhole/Spectra/090121_NM_Trypsin_20.mzXML";
		String outputFile = "/home/sangtaekim/Research/ToolDistribution/Test/HeckCIDTryp200Spec.mgf";
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));

		MzXMLSpectraIterator itr = new MzXMLSpectraIterator(specFile);
		int specCount = 0;
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			if(spec.getActivationMethod() == ActivationMethod.CID)
			{
				spec.outputMgf(out);
				specCount++;
			}
			if(specCount == 200)
				break;
		}
		out.close();
		System.out.println("Done");
	}
	
	private static void makeOMSSAScriptMerged() throws Exception
	{
		File specDir = new File("/home/sangtaekim/Research/Data/HeckWhole/MergedSpectra");
		File outputDir = new File("/home/sangtaekim/Research/Data/HeckRevision/OMSSACIDETDResults");
		String[] databases = {"/home/sangtaekim/Research/Data/HeckRevision/database/ipi.HUMAN.v3.52.target.fasta",
				"/home/sangtaekim/Research/Data/HeckRevision/database/ipi.HUMAN.v3.52.decoy.fasta"};
		
		System.out.println("#!/bin/bash");
		for(File f : specDir.listFiles())
		{
			String fileName = f.getName();
			boolean isTryp;
			if(fileName.contains("Tryp"))
				isTryp = true;
			else
				isTryp = false;
			
			if(fileName.endsWith(".mgf"))
			{
				for(int i=0; i<databases.length; i++)
				{
					String outFilePath = outputDir.getPath() + File.separator + fileName.substring(0, fileName.lastIndexOf('.')) + "_" + i + ".csv";
					System.out.println(
							"./omssacl -fm " + f.getPath() + " " +
							"-oc " + outFilePath + " " +
							"-to 0.5 -te 0.05 -tez 0 -hl 1 -he 100 -mf 3 " +
							"-d " + databases[i] + " " +
							"-ht 15 -no 7 -zcc 1 -v 2 -tem 0 -tom 0 " +
							"-i 1,2,4,5 " +
							"-e " + (isTryp ? "0" : "21"));
				}
			}
		}
	}
	
	private static void makeOMSSAScript() throws Exception
	{
		File specDir = new File("/home/sangtaekim/Research/Data/HeckWhole/MgfSpectra");
		File outputDir = new File("/home/sangtaekim/Research/Data/HeckRevision/OMSSAResults");
		String[] databases = {"/home/sangtaekim/Research/Data/HeckRevision/database/ipi.HUMAN.v3.52.target.fasta",
				"/home/sangtaekim/Research/Data/HeckRevision/database/ipi.HUMAN.v3.52.decoy.fasta"};
		
		System.out.println("#!/bin/bash");
		for(File f : specDir.listFiles())
		{
			String fileName = f.getName();
			boolean isCID;
			if(fileName.contains("CID"))
				isCID = true;
			else
				isCID = false;
			boolean isTryp;
			if(fileName.contains("Tryp"))
				isTryp = true;
			else
				isTryp = false;
			
			if(fileName.endsWith(".mgf"))
			{
				for(int i=0; i<databases.length; i++)
				{
					String outFilePath = outputDir.getPath() + File.separator + fileName.substring(0, fileName.lastIndexOf('.')) + "_" + i + ".csv";
					System.out.println(
							"./omssacl -fm " + f.getPath() + " " +
							"-oc " + outFilePath + " " +
							"-to 0.5 -te 0.05 -tez 0 -hl 1 -he 100 -mf 3 " +
							"-d " + databases[i] + " " +
							"-ht 15 -no 7 -zcc 1 -v 2 -tem 0 -tom 0 " +
							"-i " + (isCID ? "1,4" : "2,5") + " " +
							"-e " + (isTryp ? "0" : "21"));
				}
			}
		}
	}
	
	private static void mergeSpectralPairs() throws Exception
	{
		File dir = new File(System.getProperty("user.home")+"/Research/Data/HeckWhole/Spectra");
		File outputDir = new File(System.getProperty("user.home")+"/Research/Data/HeckWhole/MergedSpectra");
		for(File f : dir.listFiles())
		{
			String fileName = f.getName();
			Iterator<Spectrum> itr = null;
			if(fileName.endsWith(".mzXML"))
				itr = new MzXMLSpectraIterator(f.getPath());
			else if(fileName.endsWith(".mgf"))
				itr = new SpectraIterator(f.getPath(), new MgfSpectrumParser());
			else
				continue;
			
			System.out.print(f.getName() + ": ");
			String num = fileName.substring(fileName.lastIndexOf('_')+1, fileName.lastIndexOf('.'));
			String enzyme;
			if(fileName.contains("LysN"))
				enzyme = "LysN";
			else
				enzyme = "Tryp";
			String outputFileName = enzyme + "_" + num + "_Merged.mgf";
			PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputDir.getPath() + File.separator + outputFileName)));
			
			Spectrum prevSpec = null;
			int numPairs = 0;
			while(itr.hasNext())
			{
				Spectrum spec = itr.next();
				if(prevSpec == null)
				{
					prevSpec = spec;
					continue;
				}
				else
				{
					if((Math.round(prevSpec.getPrecursorPeak().getMass()-spec.getPrecursorPeak().getMass()) >= 0.001f) ||
							(prevSpec.getScanNum()+1 != spec.getScanNum()) 
//							|| (prevSpec.getActivationMethod() == ActivationMethod.CID && spec.getActivationMethod() == ActivationMethod.ETD 
//							|| prevSpec.getActivationMethod()== ActivationMethod.ETD && spec.getActivationMethod() == ActivationMethod.CID)
					)
					{
						prevSpec = spec;
						continue;
					}
					numPairs++;
					Spectrum mergedSpec = prevSpec.getCloneWithoutPeakList();
					mergedSpec.setActivationMethod(null);
					mergedSpec.addAll(prevSpec);
					mergedSpec.addAll(spec);
					Collections.sort(mergedSpec, new Peak.MassComparator());
					mergedSpec.outputMgf(out);
				}
			}
			System.out.println(numPairs);
			out.close();
		}
		System.out.println("Done");
	}
	
	private static HashMap<String,Float> mascotThresholds;
	private static HashMap<String,Float> msgfThresholds;
//	private static AminoAcidSet aaSet;
	static {
		mascotThresholds = new HashMap<String,Float>();
		mascotThresholds.put("CIDTryp", 39.99f);
		mascotThresholds.put("ETDTryp", 42.0f);
		mascotThresholds.put("CIDLysN", 36.75f);
		mascotThresholds.put("ETDLysN", 32.14f);
		
		msgfThresholds = new HashMap<String,Float>();
		msgfThresholds.put("CIDTryp", 2.8650214E-11f);
		msgfThresholds.put("ETDTryp", 2.7546401E-11f);
		msgfThresholds.put("CIDLysN", 7.191407E-12f);
		msgfThresholds.put("ETDLysN", 1.7512046E-11f);
//		aaSet = AminoAcidSet.getAminoAcidSet("/home/sangtaekim/Research/Data/HeckRevision/AAMassesPhospho.txt");
		
	}
	
	public static void analyzeMSGDIDs() throws Exception
	{
		String[] methods = {"CID", "ETD"};
		String[] enzymes = {"Tryp", "LysN"};
		for(String enzyme : enzymes)
			for(String method : methods)
				analyzeMSGDIDs(method, enzyme);
		
	}
	
	public static void analyzeMSGDIDs(String method, String enzymeName) throws Exception
	{
		String dir = System.getProperty("user.home")+"/Research/Data/HeckRevision/PhosMSGF/All/";
		String fileName = "MSGD_"+enzymeName+"All_0.txt";
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet("/home/sangtaekim/Research/Data/HeckRevision/AAMassesPhospho.txt");
	
		System.out.println(fileName);
		
		int numSpecs = 0;
		HashSet<String> pepSet = new HashSet<String>();
		HashSet<String> modPepSet = new HashSet<String>();
		HashSet<String> annoSet = new HashSet<String>();
		Histogram<Integer> nctHist = new Histogram<Integer>();
		Histogram<String> ptmHist = new Histogram<String>();
		
		Enzyme enzyme = Enzyme.getEnzymeByName(enzymeName);
		float threshold = msgfThresholds.get(method+enzymeName);
		BufferedLineReader in = new BufferedLineReader(dir+fileName);
		String s;
		String prevScan = "asdfasdf";
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			if(token.length != 11)
				continue;
			if(!token[2].equalsIgnoreCase(method))
				continue;
			
			if(token[1].equalsIgnoreCase(prevScan))
				continue;
			else
				prevScan = token[1];
			
			float specProb = Float.parseFloat(token[9]);
			if(specProb > threshold)
				continue;
			
			numSpecs++;
			
			String annotation = token[3];
			String pepStr = annotation.substring(annotation.indexOf('.')+1, annotation.lastIndexOf('.'));
			if(!pepSet.contains(pepStr.toUpperCase()))
			{
				pepSet.add(pepStr.toUpperCase());
				modPepSet.add(pepStr);
				annoSet.add(annotation);
			}
		}

		for(String annotation : annoSet)
		{
			int nCT = enzyme.getNumCleavedTermini(annotation, aaSet);
			nctHist.add(nCT);
		}
		
		for(String pepStr : modPepSet)
		{
			int numPhosS = 0;
			int numPhosT = 0;
			int numPhosY = 0;
			for(int i=0; i<pepStr.length(); i++)
			{
				char c = pepStr.charAt(i);
				if(c == 's')
					numPhosS++;
				else if(c == 't')
					numPhosT++;
				else if(c == 'y')
					numPhosY++;
			}
			ptmHist.add(""+numPhosS+numPhosT+numPhosY);
//			System.out.println(pepStr);
		}
		
		System.out.println("NumPeptides: " + pepSet.size());
		System.out.println("NumSpectra: " + numSpecs);
		System.out.println("Number of cleaved sites:");
		nctHist.printSorted();
		System.out.println("PTM:");
		ptmHist.printSorted();
		System.out.println();
	}
	
	public static void analyzeMascotIDs() throws Exception
	{
		String[] methods = {"CID", "ETD"};
		String[] enzymes = {"Tryp", "LysN"};
		for(String enzyme : enzymes)
			for(String method : methods)
				analyzeMascotIDs(method, enzyme);
	}
	
	public static void analyzeMascotIDs(String method, String enzymeName) throws Exception
	{
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet("/home/sangtaekim/Research/Data/HeckRevision/AAMassesPhospho.txt");
		String dir = System.getProperty("user.home")+"/Research/Data/HeckRevision/PhosMascot23/";
		String fileName = "Phos"+method+enzymeName+"Target.txt";
//		String decoyFileName = "Mascot_"+method+"_"+enzymeName+"_Decoy.txt";
	
		System.out.println(fileName);
		
		int numSpecs = 0;
		HashSet<String> pepSet = new HashSet<String>();
		HashSet<String> modPepSet = new HashSet<String>();
		HashSet<String> annoSet = new HashSet<String>();
		
		Histogram<Integer> nctHist = new Histogram<Integer>();
		Histogram<String> ptmHist = new Histogram<String>();
		
		Enzyme enzyme = Enzyme.getEnzymeByName(enzymeName);
		float threshold = mascotThresholds.get(method+enzymeName);
		BufferedLineReader in = new BufferedLineReader(dir+fileName);
		String s;
		String prevTitle = "";
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			if(token.length != 5)
				continue;
			float mascotScore = Float.parseFloat(token[3]);
			if(mascotScore <= threshold)
				continue;
			if(prevTitle.equalsIgnoreCase(token[0]))
				continue;
			else
				prevTitle = token[0];
			numSpecs++;
			
			String annotation = token[2];
			String pepStr = annotation.substring(annotation.indexOf('.')+1, annotation.lastIndexOf('.'));
			if(!pepSet.contains(pepStr.toUpperCase()))
			{
				pepSet.add(pepStr.toUpperCase());
				modPepSet.add(pepStr);
				annoSet.add(annotation);
			}
		}
		
		for(String annotation : annoSet)
		{
			int nCT = enzyme.getNumCleavedTermini(annotation, aaSet);
			nctHist.add(nCT);
		}
		
		for(String pepStr : modPepSet)
		{
			int numPhosS = 0;
			int numPhosT = 0;
			int numPhosY = 0;
			for(int i=0; i<pepStr.length(); i++)
			{
				char c = pepStr.charAt(i);
				if(c == 's')
					numPhosS++;
				else if(c == 't')
					numPhosT++;
				else if(c == 'y')
					numPhosY++;
			}
			ptmHist.add(""+numPhosS+numPhosT+numPhosY);
//			System.out.println(pepStr);
		}
		
		System.out.println("NumPeptides: " + pepSet.size());
		System.out.println("NumSpectra: " + numSpecs);
		System.out.println("Number of cleaved sites:");
		nctHist.printSorted();
		System.out.println("PTM:");
		ptmHist.printSorted();
		System.out.println();
	}
	
	public static void compareMzXMLAndMgf() throws Exception
	{
		File mgfDir = new File("/home/sangtaekim/Research/Data/HeckRevision/mgf");
		File mzXMLDir = new File("/home/sangtaekim/Research/Data/HeckRevision/spectra");

		HashMap<String,MzXMLSpectraMap> specTable = new HashMap<String,MzXMLSpectraMap>(); 
		String[] enzymes = {"Tryp", "LysN"};
		for(String enzyme : enzymes)
		{
			File dir = new File(mzXMLDir+File.separator+enzyme);
			for(File f : dir.listFiles())
			{
				if(f.getName().endsWith("mzXML"))
				{
					MzXMLSpectraMap map = new MzXMLSpectraMap(f.getPath());
					String name = f.getName().substring(0, f.getName().lastIndexOf('.'));
					name = name.toLowerCase();
					specTable.put(name, map);
				}
			}
		}
		
		int largeDiff = 0;
		for(File f : mgfDir.listFiles())
		{
			if(f.getName().endsWith("mgf"))
			{
				System.out.println(f.getName());
				SpectraIterator itr = new SpectraIterator(f.getPath(), new MgfSpectrumParser());
				while(itr.hasNext())
				{
					Spectrum spec = itr.next();
					String title = spec.getTitle();
					String[] token = title.split("\\s+");
					String fileName = token[6].substring(0, token[6].lastIndexOf('_'));
					int scanNum = Integer.parseInt(token[token.length-1]);
					float precursorMz = spec.getPrecursorPeak().getMz();
					
					Spectrum mzXMLSpec = specTable.get(fileName).getSpectrumBySpecIndex(scanNum);
					assert(mzXMLSpec != null);
					float mzXMLPrecursorMz = mzXMLSpec.getPrecursorPeak().getMz();
					
					float difference = precursorMz - mzXMLPrecursorMz;
					float diffPPM = difference/spec.getParentMass()*1e6f;
					if(diffPPM > 5)
					{	largeDiff++;
						System.out.println(difference + " " + diffPPM);
					}
				}
				System.out.println("NumLargeDiff: " + largeDiff);
			}
		}
		
	}
	
	public static void processMSGDResult() throws Exception
	{
		String[] enzymes = {"Tryp", "LysN"};
		String[] dbs = {"Target", "Decoy"};
		for(String enzyme : enzymes)
			for(String db : dbs)
				processMSGDResult(enzyme, db);
	}
	
	public static void processMSGDResult(String enzyme, String db) throws Exception
	{
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSet("/home/sangtaekim/Research/Data/HeckRevision/AAMassesPhospho.txt");
		String dir = System.getProperty("user.home")+"/Research/Data/HeckRevision/MSGPSearch/";
		String fileName = "MSGD_"+enzyme+"_"+db+"_All.txt";
		String outputFileName = "MSGD_"+enzyme+"_"+db+".txt";
		
		BufferedLineReader in = new BufferedLineReader(dir+"/"+fileName);
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(dir+outputFileName)));

//		HashMap<String,MzXMLSpectraMap> specMap = new HashMap<String,MzXMLSpectraMap>();
		
		String s;
		int prevScanNum = -1;
		float bestSpecProb = 1;
		String bestStr = null;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#"))
			{
				out.println(s);
				continue;
			}
			String[] token = s.split("\t");
			if(token.length != 11)
				continue;
			
			String specFile = token[0].split(":")[0];
			int scanNum = Integer.parseInt(token[1]);
			
//			MzXMLSpectraMap map = specMap.get(specFile);
//			if(map == null)
//			{
//				map = new MzXMLSpectraMap(specFile);
//				specMap.put(specFile, map);
//			}
//			Spectrum spec = map.getSpectrumByScanNum(scanNum);
//			assert(spec != null);
//			float parentMass = spec.getParentMass();
			
			String pep = token[3];
			pep = pep.substring(pep.indexOf('.')+1, pep.lastIndexOf('.'));
			Peptide peptide = aaSet.getPeptide(pep);
			float theoMass = peptide.getParentMass();
			float massError = Float.parseFloat(token[10]);
			float massErrorPPM = massError/theoMass*1e6f;
//			if(massErrorPPM > 5)
//				continue;
			
			float specProb = Float.parseFloat(token[9]);
			if(scanNum == prevScanNum)
			{
				if(specProb < bestSpecProb)
				{
					bestSpecProb = specProb;
					bestStr = s;
				}
			}
			else	// new scan
			{
				if(bestStr != null)
					out.println(bestStr);
				
				prevScanNum = scanNum;
				bestSpecProb = specProb;
				bestStr = s;
			}
		}
		if(bestStr != null)
			out.println(bestStr);
		
		in.close();
		out.close();
		System.out.println("Done");
	}
	
	public static void splitMascotResult() throws Exception
	{
		String[] methods = {"CID", "ETD"};
		String[] enzymes = {"Tryp", "LysN"};
		for(String method : methods)
			for(String enzyme : enzymes)
				splitMascotResult(method, enzyme);
	}
	
	public static void splitMascotResult(String method, String enzyme) throws Exception
	{
		String dir = System.getProperty("user.home")+"/Research/Data/HeckRevision/MascotResults/";
		String fileName = "Mascot_"+method+"_"+enzyme+".txt";
		String targetFileName = "Mascot_"+method+"_"+enzyme+"_Target.txt";
		String decoyFileName = "Mascot_"+method+"_"+enzyme+"_Decoy.txt";
		
		PrintStream targetOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(dir+targetFileName)));
		PrintStream decoyOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(dir+decoyFileName)));
		PrintStream out = null;
		
		HashSet<String> targetTitles = new HashSet<String>();
		HashSet<String> decoyTitles = new HashSet<String>();
		
		BufferedLineReader in = new BufferedLineReader(dir+"/"+fileName);
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			if(token.length != 5)
				continue;
			if(token[4].contains("_reversed"))
			{
				out = decoyOut;
				if(decoyTitles.contains(token[0]))
					continue;
				decoyTitles.add(token[0]);
			}
			else
			{
				out = targetOut;
				if(targetTitles.contains(token[0]))
					continue;
				targetTitles.add(token[0]);
			}
			out.println(s);
		}
		targetOut.close();
		decoyOut.close();
		
		int targetOnly = 0, decoyOnly = 0, shared = 0;
		for(String title : targetTitles)
		{
			if(decoyTitles.contains(title))
			{
				shared++;
				decoyTitles.remove(title);
			}
			else
				targetOnly++;
		}
		decoyOnly = decoyTitles.size();
		System.out.println("TargetOnly: " + targetOnly);
		System.out.println("DecoyOnly: " + decoyOnly);
		System.out.println("Shared: " + shared);
		System.out.println("Done");
	}
	
	public static void splitTargetDecoy() throws Exception
	{
		String dir = System.getProperty("user.home")+"/Research/Data/HeckRevision/database/";
		String fileName = "ipi.HUMAN.v3.52.smartrev.fasta";
		String targetFileName = "ipi.HUMAN.v3.52.target.fasta";
		String decoyFileName = "ipi.HUMAN.v3.52.decoy.fasta";
		
		PrintStream targetOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(dir+targetFileName)));
		PrintStream decoyOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(dir+decoyFileName)));
		PrintStream out = null;
		
		BufferedLineReader in = new BufferedLineReader(dir+"/"+fileName);
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith(">"))
			{
//				if(s.contains("_reversed"))
				if(s.contains("IPI:REV"))
					out = decoyOut;
				else
					out = targetOut;
			}
			if(out != null)
				out.println(s);
		}
		
		targetOut.close();
		decoyOut.close();
		System.out.println("Done");
	}
	
	public static void analyzeMascotResults() throws Exception
	{
		String fileName = "/home/sangtaekim/Research/Data/HeckRevision/MascotResults/Mascot_CID_Tryp.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		
		int numSpectra = 0;
		int numS = 0;
		int numT = 0;
		int numY = 0;
		int numDouble = 0;
		HashSet<String> pepSet = new HashSet<String>();
		int[] numPTMHist = new int[100];
		in.readLine();
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length != 5)
				continue;
			
			String title = token[0];
			int charge = Integer.parseInt(token[1]);
			String annotation = token[2];
			float mascotScore = Float.parseFloat(token[3]);
			if(mascotScore < 30)
				continue;
			String peptide = annotation.substring(annotation.indexOf('.')+1, annotation.lastIndexOf('.'));
			if(peptide.contains("B") || peptide.contains("J") || peptide.contains("O")
					|| peptide.contains("U") || peptide.contains("X") || peptide.contains("Z"))
				continue;
			pepSet.add(peptide.toUpperCase());
			numSpectra++;
			if(peptide.contains("s"))
				numS++;
			if(peptide.contains("t"))
				numT++;
			if(peptide.contains("y"))
				numY++;
			int numPTMs = 0;
			for(int i=0; i<peptide.length(); i++)
			{
				char aa = peptide.charAt(i);
				if(aa == 's' || aa == 't' || aa == 'y')
					numPTMs++;
			}
			numPTMHist[numPTMs]++;
			String protein = token[4];
		}
		
		System.out.println("NumSpectra: " + numSpectra);
		System.out.println("NumDistinctPeptides: " + pepSet.size());
		System.out.println("NumPhosS: " + numS);
		System.out.println("NumPhosT: " + numT);
		System.out.println("NumPhosY: " + numY);
		System.out.println("NumMultiplePTMs: " + numDouble);
		System.out.println("NumPTMsHist:");
		for(int i=0; i<numPTMHist.length; i++)
		{
			if(numPTMHist[i] > 0)
				System.out.println(i+"\t"+numPTMHist[i]);
		}
	}
}
