package misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;

import parser.BufferedLineReader;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraIterator;
import msgf.AminoAcidGraph;
import msgf.GeneratingFunction;
import msgf.ProfileGF;
import msgf.ScoredSpectrum;
import msgf.Tolerance;
import msscorer.NewScorerFactory;
import msutil.ActivationMethod;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.Enzyme;
import msutil.IonType;
import msutil.Peptide;
import msutil.ScoredString;
import msutil.Sequence;
import msutil.SpectraContainer;
import msutil.SpectraIterator;
import msutil.Spectrum;
import msutil.SpectrumAnnotator;
import msutil.WindowFilter;

public class CIDETDPairs {
	public static void main(String argv[]) throws Exception
	{
//		mergeGRIDResults();
//		rescore();
//		intersection();
//		better();
		vennDiagram();
//		makeAnnotatedSpectra();
//		computeCoverage();
//		computeCoveragePairs();
//		testGappedPeptidePerformancePair(0.96f, 0.4f);
//		processIonProb();
//		generateRankDist();
//		generateErrorDist();
//		sortResult(new File("/home/sangtaekim/Research/MSGF2D/rescored/old"));
//		addScanNumDir(new File("/home/sangtaekim/Research/Data/HeckWhole/MgfSpectra"));
//		splitConcatenatedResults(new File("/home/sangtaekim/Research/Data/HeckWhole/ResultsForNuno"), true);
//		splitConcatenatedResults(new File("/home/sangtaekim/Research/Data/ISBETD/MSGFDB0316"), true); 
//		countPairs("/home/sangtaekim/Research/Data/HeckWhole/Spectra");
	}
	
	
	public static void countPairs(String dirName) throws Exception
	{
		int numPairs = 0;
		int[][] numPairsC = new int[2][3];
		File dir = new File(dirName);
		for(File f : dir.listFiles())
		{
			int enzIndex = 0;
			if(f.getName().contains("LysN"))
				enzIndex = 1;
			if(f.getName().endsWith(".mzXML"))
			{
				MzXMLSpectraIterator itr = new MzXMLSpectraIterator(f.getPath());
				float prevPrecursorMz = 0;
				int prevCharge = 0;
				while(itr.hasNext())
				{
					Spectrum spec = itr.next();
					float precursorMz = spec.getPrecursorPeak().getMz();
					int charge = spec.getCharge();
					if(precursorMz == prevPrecursorMz && charge == prevCharge)
					{
						numPairs++;
						if(charge == 2)
							numPairsC[enzIndex][0]++;
						else if(charge == 3)
							numPairsC[enzIndex][1]++;
						else if(charge > 3)
							numPairsC[enzIndex][2]++;
					}
					prevPrecursorMz = precursorMz;
					prevCharge = charge;
				}
			}
			else if(f.getName().endsWith(".mgf"))
			{
				SpectraIterator itr = new SpectraIterator(f.getPath(), new MgfSpectrumParser());
				float prevPrecursorMz = 0;
				int prevCharge = 0;
				while(itr.hasNext())
				{
					Spectrum spec = itr.next();
					float precursorMz = spec.getPrecursorPeak().getMz();
					int charge = spec.getCharge();
					if(precursorMz == prevPrecursorMz && charge == prevCharge)
					{
						numPairs++;
						if(charge == 2)
							numPairsC[enzIndex][0]++;
						else if(charge == 3)
							numPairsC[enzIndex][1]++;
						else if(charge > 3)
							numPairsC[enzIndex][2]++;
					}
					prevPrecursorMz = precursorMz;
					prevCharge = charge;
				}
			}
		}
		System.out.println("NumPairs: " + numPairs);
		System.out.println("NumPairs Tryp C2: " + numPairsC[0][0]);
		System.out.println("NumPairs Tryp C3: " + numPairsC[0][1]);
		System.out.println("NumPairs Tryp C4+: " + numPairsC[0][2]);
		System.out.println("NumPairs LysN C2: " + numPairsC[1][0]);
		System.out.println("NumPairs LysN C3: " + numPairsC[1][1]);
		System.out.println("NumPairs LysN C4+: " + numPairsC[1][2]);
	}
	
	public static void splitConcatenatedResults(File dir, boolean bestOnly) throws Exception
	{
		if(!dir.isDirectory())
			return;
		for(File f : dir.listFiles())
		{
//			if(f.getName().endsWith("_0.txt"))
			if(f.getName().endsWith("_Target.txt"))
			{
				/////////
				String name = f.getName();
				name = name.substring(name.indexOf('_')+1, name.indexOf('_', name.indexOf('_')+1));
				String num = f.getName().substring(f.getName().lastIndexOf("Target.txt")-3, f.getName().lastIndexOf("_Target.txt"));
				String targetName = name+"_"+num+".txt";
				//////////////
//				String name = f.getName().substring(0, f.getName().lastIndexOf("_0.txt"));
//				String targetName = name + "_Target.txt";
				PrintStream targetOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(f.getParent()+File.separator+targetName)));
				String decoyName = name + "_Decoy.txt";
				PrintStream decoyOut = new PrintStream(new BufferedOutputStream(new FileOutputStream(f.getParent()+File.separator+decoyName)));
				BufferedLineReader in = new BufferedLineReader(f.getPath());
				String s;
				String prevScanNum = "asdfsdaf";
				targetOut.println("FileName\tScanNum\tActivationMethod\tPrecursorMz\tCharge\tPeptide\tProtein\tOptimalScore\tPeptideScore\tSpecProb");
				while((s=in.readLine()) != null)
				{
					if(s.startsWith("#"))
						continue;
					String[] token = s.split("\t");
					if(token.length < 7)
						continue;
//					if(token[1].contains("-"))
//						continue;
					String scanNum = token[1];
					if(bestOnly && scanNum.equalsIgnoreCase(prevScanNum))
						continue;
					else
						prevScanNum = scanNum;
					if(token[6].contains("IPI:REV") || token[6].contains("DECOY"))
						decoyOut.println(s);
					else
						targetOut.println(s);
				}
				in.close();
				targetOut.close();
				decoyOut.close();
			}
		}
		System.out.println("Done");
	}
	
	public static void addScanNumDir(File dir) throws Exception
	{
		if(!dir.isDirectory())
			return;
		for(File f : dir.listFiles())
		{
			if(!f.getName().endsWith(".mgf"))
				continue;
			addScanNum(f);
		}
		System.out.println("Done");
	}
	
	public static void addScanNum(File f) throws Exception
	{
		String output = f.getName().substring(f.getName().indexOf("NM_")+3);
		output = f.getParent()+File.separator+output;
		System.out.println(output);

		SpectraIterator itr = new SpectraIterator(f.getPath(), new MgfSpectrumParser());
		SpectraContainer container = new SpectraContainer();
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			String title = spec.getTitle();
			int index = title.indexOf("ScanNumber: ")+"ScanNumber: ".length();
			String scanStr = title.substring(index);
			int scanNum = Integer.parseInt(scanStr);
			spec.setScanNum(scanNum);
			container.add(spec);
		}
		container.outputMgfFile(output);
	}
	
	public static void sortResult(File dir) throws Exception
	{
		for(File f : dir.listFiles())
		{
			System.out.println(f.getName());
			if(f.getName().endsWith(".txt"))
			{
				sortResult(f.getPath(), "/home/sangtaekim/Research/MSGF2D/rescored/"+f.getName());
			}
		}
	}
	
	public static void sortResult(String fileName, String outputFileName) throws Exception
	{
		PrintStream out = null;
		if(outputFileName == null)
			out = System.out;
		else
		{
			try {
				out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFileName)));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		boolean labelExists = false;
		ArrayList<ScoredString> results = new ArrayList<ScoredString>();
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("SpecNum"))
			{
				out.println("#"+s);
				labelExists = true;
				continue;
			}
			else if(s.startsWith("#"))
				continue;
			else
			{
				if(!labelExists)
				{
					out.println("#SpecNum\tPrecursorMz\tCharge\tMSGFScore\tPeptide\tPeptideScore\tSpecProb");
					labelExists = true;
				}
				String[] token = s.split("\t");
				if(token.length != 7)
					continue;
				results.add(new ScoredString(s, Float.parseFloat(token[0])));
			}
		}
		Collections.sort(results);
		for(ScoredString ss : results)
			out.println(ss.getStr());
		out.close();
	}
	
	public static void generateErrorDist() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/MSGF2D/annotatedSpectra/ionstat_CID_Tryp_filtered.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		String[] label = null;
		int[][] table = null; 
		int[] noise = null;
		int maxRank = 0;
		int numSpecs = 0;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#MaxRanks"))
				maxRank = Integer.parseInt(s.split("\t")[1]);
			else if(s.startsWith("#RankDist"))
				numSpecs = Integer.parseInt(s.split("\t")[4]);
			String[] token = s.split("\t");
			if(token.length < 3)
				continue;
			if(token[0].equalsIgnoreCase("Rank"))
			{
				label = Arrays.copyOfRange(token, 2, token.length);
				table = new int[maxRank][label.length];
				noise = new int[maxRank];
			}
			if(token[0].length() > 0 && Character.isDigit(token[0].charAt(0)) && token[1].equalsIgnoreCase("sum"))
			{
				int rank = Integer.parseInt(token[0]);
				for(int i=2; i<token.length; i++)
				{
					if(!label[i-2].startsWith("noise"))
						table[rank-1][i-2] += Integer.parseInt(token[i]);
					else
						noise[rank-1]+=Integer.parseInt(token[i]);
				}
			}
		}
		int numNoiseBins = 0;
		for(int i=0; i<label.length; i++)
		{
			label[i] = label[i].replace("-", "Minus");
			label[i] = label[i].replace("+", "Plus");
			if(!label[i].startsWith("noise"))
			{
				System.out.print(label[i]+"=[");
				for(int j=0; j<maxRank-1; j++)
					System.out.print(table[j][i]/(float)numSpecs+" ");
				System.out.println(table[maxRank-1][i]/(float)numSpecs+"];");
			}
			else
				numNoiseBins++;
		}
		System.out.print("noise=[");
		for(int j=0; j<maxRank-1; j++)
			System.out.print(noise[j]/(float)numNoiseBins/(float)numSpecs+" ");
		System.out.println(noise[maxRank-1]/(float)numNoiseBins/(float)numSpecs+"];");
		
		System.out.println("hold off;");
		String[] color = {"b", "g", "r", "c", "m", "b*", "g*", "r*", "c*", "m*"};
		boolean holdOn = false;
		for(int i=0; i<label.length; i++)
		{
			if(!label[i].startsWith("noise"))
			{
				System.out.println("plot("+label[i]+",'"+color[i]+"','DisplayName','"+label[i].replace("Plus", "+").replace("Minus", "-")+"');");
				if(!holdOn)
				{
					System.out.println("hold on");
					holdOn = true;
				}
						
			}
		}
		System.out.println("plot("+"noise"+",'"+"k"+"','DisplayName','"+"noise"+"');");
		System.out.println("xlim([1 100]);");
		System.out.println("xlabel('Rank','FontSize',12);");
		System.out.println("ylabel('Probability','FontSize',12);");		
	}
	
	public static void generateRankDist() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/HeckWhole/AnnotatedSpectra/configForIonstat/ionstat_CID_Tryp_c2.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		String[] label = null;
		int[][] table = null; 
		int[] noise = null;
		int[] sumExplained = null;
		ArrayList<ScoredString> ionOccSum = new ArrayList<ScoredString>();
		
		String title = null;
		String method = null;
		if(fileName.contains("CID"))
			method = "CID";
		else if(fileName.contains("ETD"))
			method = "ETD";
		String enzyme = null;
		if(fileName.contains("Tryp"))
			enzyme = "Tryp";
		else if(fileName.contains("LysN"))
			enzyme = "LysN";
		String charge = null;
		if(fileName.contains("c2"))
			charge = "2";
		else if(fileName.contains("c3"))
			charge = "3";
		title = method+"-"+enzyme+" charge " + charge;
		
		int maxRank = 0;
		int numSpecs = 0;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#MaxRanks"))
				maxRank = Integer.parseInt(s.split("\t")[1]);
			else if(s.startsWith("#RankDist"))
				numSpecs = Integer.parseInt(s.split("\t")[4]);
			String[] token = s.split("\t");
			if(token.length < 3)
				continue;
			if(token[0].equalsIgnoreCase("Rank"))
			{
				label = Arrays.copyOfRange(token, 2, token.length);
				table = new int[maxRank][label.length];
				noise = new int[maxRank];
				sumExplained = new int[maxRank];
			}
			if(token[0].length() > 0 && Character.isDigit(token[0].charAt(0)) && token[1].equalsIgnoreCase("sum"))
			{
				int rank = Integer.parseInt(token[0]);
				for(int i=2; i<token.length; i++)
				{
					if(!label[i-2].startsWith("noise"))
					{
						table[rank-1][i-2] += Integer.parseInt(token[i]);
						sumExplained[rank-1] += Integer.parseInt(token[i]);
					}
					else
						noise[rank-1]+=Integer.parseInt(token[i]);
				}
			}
		}
		
		int numNoiseBins = 0;
		for(int i=0; i<label.length; i++)
		{
			label[i] = label[i].replace("-", "Minus");
			label[i] = label[i].replace("+", "Plus");
			if(!label[i].startsWith("noise"))
			{
				System.out.print(label[i]+"=[");
				float sum50 = 0;
				for(int j=0; j<maxRank-1; j++)
				{
					System.out.print(table[j][i]/(float)numSpecs+" ");
//					if(j <= 50)
						sum50 += table[j][i];
				}
				System.out.println(table[maxRank-1][i]/(float)numSpecs+"];");
				ionOccSum.add(new ScoredString(label[i], sum50));
			}
			else
				numNoiseBins++;
		}
		Collections.sort(ionOccSum, Collections.reverseOrder());
		System.out.print("noise=[");
		for(int j=0; j<maxRank-1; j++)
			System.out.print(noise[j]/(float)numNoiseBins/(float)numSpecs+" ");
		System.out.println(noise[maxRank-1]/(float)numNoiseBins/(float)numSpecs+"];");

		System.out.print("unexplained=[");
		for(int j=0; j<maxRank-1; j++)
			System.out.print(1-sumExplained[j]/(float)numSpecs+" ");
		System.out.println(1-sumExplained[maxRank-1]/(float)numSpecs+"];");
		
		System.out.println("hold off;");
		String[] color = {"b", "g", "r", "c", "m", "b*", "g*", "r*", "c*", "m*"};
		boolean holdOn = false;
		
		for(int col = 0; col<color.length; col++)
		{
			for(int i=0; i<label.length; i++)
			{
//				if(!label[i].startsWith("noise"))
				if(label[i] == ionOccSum.get(col).getStr())
				{
					System.out.println("%"+label[i]+" "+(ionOccSum));
					System.out.println("plot("+label[i]+",'"+color[col]+"','DisplayName','"+label[i].replace("Plus", "+").replace("Minus", "-")+"');");
					if(!holdOn)
					{
						System.out.println("hold on");
						holdOn = true;
					}
				}
			}
		}
//		System.out.println("plot("+"noise"+",'"+"k"+"','DisplayName','"+"noise"+"');");
		System.out.println("plot("+"unexplained"+",'"+"k"+"','DisplayName','"+"unexplained"+"');");
		System.out.println("xlim([1 100]);");
		System.out.println("ylim([0 1]);");
		System.out.println("xlabel('Rank','FontSize',14);");
		System.out.println("ylabel('Probability','FontSize',14);");
		System.out.println("title('"+title+"','FontSize',16);");
	}

	public static void processIonProb() throws Exception
	{
		String dirName = System.getProperty("user.home")+"/Research/Data/HeckWhole/AnnotatedSpectra/configForIonstat";
		File dir = new File(dirName);
		for(File f : dir.listFiles())
		{
			if(f.getName().startsWith("prob"))
			{
				processIonProb(f.getPath());
			}
		}
	}
	
	public static void processIonProb(String fileName) throws Exception
	{
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		int mode = 0;	// 1: ion, 2: noise
		
		float probUnexplained = 1;
		int numNoise = 0;
		float noiseProbSum = 0;
		ArrayList<Float> probList = new ArrayList<Float>();
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("b-H\t") || s.startsWith("z-H\t"))
				continue;
			if(mode==0)
			{
				if(s.startsWith("Ion"))
				{
					mode = 1;
					continue;
				}
			}
			else if(mode==1 && s.startsWith("noise"))
				mode = 2;
			
			if(mode == 1)
			{
				probList.add(Float.parseFloat(s.split("\t")[1]));
			}
			else if(mode == 2)
			{
				noiseProbSum += Float.parseFloat(s.split("\t")[1]);
				numNoise++;
			}
		}
		probList.add(noiseProbSum/numNoise);
		System.out.print(fileName.substring(fileName.lastIndexOf('/')+1,fileName.lastIndexOf('.'))+"=");
		System.out.print("[");
		for(int i=0; i<probList.size()-1; i++)
		{
			System.out.print(probList.get(i)+" ");
			probUnexplained -= probList.get(i);
		}
		System.out.println(probList.get(probList.size()-1)+"];");
	}
	
    public static void testGappedPeptidePerformancePair(float templateFraction, float minProbThreshold)
	{
//    	boolean printIndividualSpecResults = true;
//		System.out.println("**************** "+templateFraction+"\t"+minProbThreshold);
//		float specProbThreshold = (float)0.1/146166984;
//		int MIN_GAPPED_PEPTIDE_LENGTH = 5;
//		
//		String specFileCID = System.getProperty("user.home")+"/Research/MSGF2D/annotatedSpectra/annotated_CID_Tryp_pair.mgf";
//		String specFileETD = System.getProperty("user.home")+"/Research/MSGF2D/annotatedSpectra/annotated_ETD_Tryp_pair.mgf";
//		SpectraIterator iteratorCID = null;
//		SpectraIterator iteratorETD = null;
//		try {
//			iteratorCID = new SpectraIterator(specFileCID, new MgfSpectrumParser());
//			iteratorETD = new SpectraIterator(specFileETD, new MgfSpectrumParser());
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
//		
//		int specNum = 0;
//		int numTrueDeNovo = 0;
//		int sumDeNovoRecs = 0;
//		int[] numProcessedSpec = new int[100];
//		int[] numSpecWithCorrectGP = new int[100];
//		int[] sumSizeGP = new int[100];
//		int[] numSpectraWithLongGP = new int[100];
//		
//		int[] specProbDist = new int[100];
//		
//		if(printIndividualSpecResults)
//			System.out.println("SpecNum\tPeptide\tMSGFScore\tPeptideScore\tSpecProb\tGPLength\tIsGPCorrect");
//		while(iteratorCID.hasNext())
//		{
//			specNum++;
//			Spectrum specCID = iteratorCID.next();
//			Spectrum specETD = iteratorETD.next();
//
//			specCID.correctParentMass();
//			specETD.correctParentMass();
//			if(specCID.getCharge() != 2)
//				continue;
//
//			Peptide annotation = specCID.getAnnotation();
//			int length = annotation.size();
//			
//			AminoAcidGraph graph = new AminoAcidGraph(specCID.getParentMass(), AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys());
//
//			if(graph.getPMNode() == null)
//				continue;
//
//			ScoredSpectrum<NominalMass> scoredSpecCID = NewScorerFactory.get(ActivationMethod.CID, Enzyme.TRYPSIN).getScoredSpectrum(specCID);
//			ScoredSpectrum<NominalMass> scoredSpecETD = NewScorerFactory.get(ActivationMethod.CID, Enzyme.TRYPSIN).getScoredSpectrum(specCID);
//			
//			GeneratingFunction<NominalMass> gf = new GeneratingFunction<NominalMass>(new ScoredSpectrumSumPairs<NominalMass>(scoredSpecCID, scoredSpecETD), graph);
//			graph.allowNonEnzymaticCleavage();
//			gf.computeGeneratingFunction();
//			
//			ProfileGF<NominalMass> profGf = new ProfileGF<NominalMass>(gf);
//			profGf.computeProfile(specProbThreshold);
//			Sequence<NominalMass> gappedPeptide = profGf.getGappedPeptideWithNominalMasses(templateFraction, minProbThreshold);
//			
//			float specProb = gf.getSpectralProbability(annotation);
//			int specProbIndex = Math.round((float)(-Math.log10(specProb)));
//			specProbDist[specProbIndex]++;
//			boolean isGPTrue = gappedPeptide.isMatchedToNominalMasses(annotation, false);
//			if(gappedPeptide.size() >= MIN_GAPPED_PEPTIDE_LENGTH)
//			{
//				numSpectraWithLongGP[length]++;		
//				if(isGPTrue)
//					numSpecWithCorrectGP[length]++;
//				sumSizeGP[length] += gappedPeptide.size();
//			}
//			else
//				isGPTrue = false;
//			int deNovoScore = gf.getMaxScore()-1;
//			int peptideScore = gf.getScore(annotation);
//			numProcessedSpec[length]++;
//			if(deNovoScore == peptideScore)
//				numTrueDeNovo++;
//			int numDeNovoSequences = (int) gf.getNumEqualOrBetterPeptides(deNovoScore);
//			assert(specProb > 0);
//			sumDeNovoRecs += numDeNovoSequences;
//			if(printIndividualSpecResults)
//				System.out.println(specNum+"\t"+annotation+"\t"+deNovoScore+"\t"+peptideScore+"\t"+specProb+"\t"+gappedPeptide.size()+"\t"+isGPTrue+"\t"+numDeNovoSequences);
//		}
//		
//		int overallNumProcessedSpec = 0;
//		int overallNumSpectraWithLongGP = 0;
//		int overallNumSpecWithCorrectGP = 0;
//		int overallSumSizeGP = 0;
//		
//		for(int i=0; i<numProcessedSpec.length; i++)
//		{
//			if(numProcessedSpec[i] > 0)
//			{
//				overallNumProcessedSpec += numProcessedSpec[i];
//				overallNumSpectraWithLongGP += numSpectraWithLongGP[i];
//				overallNumSpecWithCorrectGP += numSpecWithCorrectGP[i];
//				overallSumSizeGP += sumSizeGP[i];
//				System.out.println(i+"\t"+numProcessedSpec[i]+"\t"+
//						(numSpectraWithLongGP[i]/(float)numProcessedSpec[i])+"\t"+
//						(numSpecWithCorrectGP[i]/(float)numSpectraWithLongGP[i])+"\t"+
//						(sumSizeGP[i]/(float)numSpectraWithLongGP[i]));
//			}
//		}
//		
//		System.out.println("Overall");
//		System.out.println("NumSpectra\t"+overallNumProcessedSpec);
//		System.out.println("Sensitivity(#SpecWithLongGP/#Spectra)\t"+(overallNumSpectraWithLongGP/(float)overallNumProcessedSpec));
//		System.out.println("Accuracy(#CorrectLongGP/#SpecWithLongGP)\t"+(overallNumSpecWithCorrectGP/(float)overallNumSpectraWithLongGP));
//		System.out.println("AverageLongGappedPeptideLength\t"+(overallSumSizeGP/(float)overallNumSpectraWithLongGP));
//		System.out.println();
	}
    
	public static void computeCoveragePairs() throws Exception
	{
		String fileNameCID = System.getProperty("user.home")+"/Research/MSGF2D/annotatedSpectra/annotated_CID_LysN_pair.mgf";
		String fileNameETD = System.getProperty("user.home")+"/Research/MSGF2D/annotatedSpectra/annotated_ETD_LysN_pair.mgf";
		Tolerance tolerance = new Tolerance(0.2f, false);
		WindowFilter filter = new WindowFilter(6, 50);
		SpectraIterator itrCID = new SpectraIterator(fileNameCID, new MgfSpectrumParser());
		SpectraIterator itrETD = new SpectraIterator(fileNameETD, new MgfSpectrumParser());
		
		IonType[] ionTypesCID = {IonType.Y, IonType.B, IonType.A};
		IonType[] ionTypesETD = {IonType.getIonType("z+H"), IonType.C, IonType.Y, IonType.B};
		
		int numSpecs = 0;
		float covSum = 0;
		String id;
		if(fileNameCID.contains("Tryp"))
			id = "Tryp";
		else
			id = "LysN";
		System.out.print(id+"=[");
		while(itrCID.hasNext())
		{
			Spectrum specCID = itrCID.next();
			Spectrum specETD = itrETD.next();
			numSpecs++;
			specCID = filter.apply(specCID);
			specETD = filter.apply(specETD);
			Peptide annotation = specCID.getAnnotation();
			SpectrumAnnotator annotatorCID = new SpectrumAnnotator(specCID, annotation);
			SpectrumAnnotator annotatorETD = new SpectrumAnnotator(specETD, annotation);
			
			HashSet<Integer> indices = new HashSet<Integer>();
			for(Integer index : annotatorCID.getCoveredCleavages(ionTypesCID, tolerance))
				indices.add(index);
			for(Integer index : annotatorETD.getCoveredCleavages(ionTypesETD, tolerance))
				indices.add(index);
			float coverage = indices.size()/(float)(annotation.size()-1);
			System.out.print(coverage+" ");
			covSum += coverage;
		}
		System.out.println("];");
		System.out.println(fileNameCID);
		System.out.print("Ion (CID): ");
		for(IonType ion : ionTypesCID)
			System.out.print("\t"+ion.getName());
		System.out.println();
		System.out.println(fileNameETD);
		System.out.print("Ion (ETD): ");
		for(IonType ion : ionTypesETD)
			System.out.print("\t"+ion.getName());
		System.out.println();
		System.out.println("NumSpecs: " + numSpecs);
		System.out.println("Coverage: " + covSum/numSpecs);
	}
	
	public static void computeCoverage() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/MSGF2D/annotatedSpectra/annotated_ETD_LysN.mgf";
		Tolerance tolerance = new Tolerance(0.5f, false);
		WindowFilter filter = new WindowFilter(6, 50);
		SpectraIterator itr = new SpectraIterator(fileName, new MgfSpectrumParser());
//		IonType[] ionTypes = {IonType.Y, IonType.B, IonType.A};
		IonType[] ionTypes = {IonType.getIonType("z+H"), IonType.C, IonType.Y, IonType.B};
		int numSpecs = 0;
		float covSum = 0;
		String id = fileName.substring(fileName.indexOf('_')+1, fileName.lastIndexOf('.'));
		System.out.print(id+"=[");
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			numSpecs++;
			spec = filter.apply(spec);
			Peptide annotation = spec.getAnnotation();
			SpectrumAnnotator annotator = new SpectrumAnnotator(spec, annotation);
			float coverage = annotator.getCoverage(ionTypes, tolerance);
//			System.out.println(spec.getScanNum()+"\t"+spec.getPrecursorPeak().getMz()+"\t"+spec.getCharge()+
//					"\t"+annotation+"\t"+coverage);
			System.out.print(coverage+" ");
			covSum += coverage;
		}
		System.out.println("];");
		System.out.println(fileName);
		System.out.print("Ion: ");
		for(IonType ion : ionTypes)
			System.out.print("\t"+ion.getName());
		System.out.println();
		System.out.println("NumSpecs: " + numSpecs);
		System.out.println("Coverage: " + covSum/numSpecs);
	}
	
	public static void makeAnnotatedPairedSpectra() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/MSGF2D/rescored/LysNSum_0.txt";
		
		String specFileCID1 = System.getProperty("user.home")+"/Research/MSGF2D/CID_paired_01_LysN.mgf";
		String specFileCID2 = System.getProperty("user.home")+"/Research/MSGF2D/CID_paired_03_LysN.mgf";
		String specFileETD1 = System.getProperty("user.home")+"/Research/MSGF2D/ETD_paired_01_LysN.mgf";
		String specFileETD2 = System.getProperty("user.home")+"/Research/MSGF2D/ETD_paired_03_LysN.mgf";
		
		String outputFileCID = System.getProperty("user.home")+"/Research/MSGF2D/annotatedSpectra/annotated_CID_LysN_pair.mgf";
		String outputFileETD = System.getProperty("user.home")+"/Research/MSGF2D/annotatedSpectra/annotated_ETD_LysN_pair.mgf";
		
		SpectraContainer newContainerCID = new SpectraContainer();
		SpectraContainer newContainerETD = new SpectraContainer();
//		float threshold = 1.9050785E-11f;	// Tryp
		float threshold = 1.8293823E-11f;	// LysN
		
		String s;
		Hashtable<String, String> pepTable = new Hashtable<String, String>();
		BufferedLineReader in = new BufferedLineReader(fileName);
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("SpecNum") || s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			if(token.length != 7)
				continue;
			float specProb = Float.parseFloat(token[6]);
			if(specProb < threshold)
			{
				String pep = token[4].substring(token[4].indexOf('.')+1, token[4].lastIndexOf('.'));
				if(pepTable.get(pep) == null)
					pepTable.put(pep, s);
				else
				{
					String[] token2 = pepTable.get(pep).split("\t");
					assert(token2.length == 7);
					float existingSpecProb = Float.parseFloat(token2[6]);
					if(specProb < existingSpecProb)
						pepTable.put(pep, s);
				}
			}
		}
		
		SpectraContainer containerCID1 = new SpectraContainer(specFileCID1, new MgfSpectrumParser());
		SpectraContainer containerCID2 = new SpectraContainer(specFileCID2, new MgfSpectrumParser());
		SpectraContainer containerETD1 = new SpectraContainer(specFileETD1, new MgfSpectrumParser());
		SpectraContainer containerETD2 = new SpectraContainer(specFileETD2, new MgfSpectrumParser());
		
		for(String str : pepTable.values())
		{
			String[] token = str.split("\t");
			int specNum = Integer.parseInt(token[0]);
//			float precursorMz = Float.parseFloat(token[1]);
			String annotation = token[4].substring(token[4].indexOf('.')+1, token[4].lastIndexOf('.'));
			int charge = Integer.parseInt(token[2]);
			float precursorMz = (new Peptide(annotation).getParentMass()+charge*(float)Composition.PROTON)/charge;
			Spectrum specCID1 = null;
			Spectrum specCID2 = null;
			Spectrum specETD1 = null;
			Spectrum specETD2 = null;
			if(specNum-1 < containerCID1.size())
			{
				specCID1 = containerCID1.get(specNum-1);
				specETD1 = containerETD1.get(specNum-1);
			}
			if(specNum-1 < containerCID2.size())
			{
				specCID2 = containerCID2.get(specNum-1);
				specETD2 = containerETD2.get(specNum-1);
			}
			
			Spectrum specCID = null;
			Spectrum specETD = null;
			if(specCID1 == null)
			{
				specCID = specCID2;
				specETD = specETD2;
			}
			else if(specCID2 == null)
			{
				specCID = specCID1;
				specETD = specETD1;
			}
			else
			{
				if(Math.abs(precursorMz-specCID1.getPrecursorPeak().getMz()) < Math.abs(precursorMz-specCID2.getPrecursorPeak().getMz()))
				{
					specCID = specCID1;
					specETD = specETD1;
				}
				else
				{
					specCID = specCID2;
					specETD = specETD2;
				}
				if(Math.abs(precursorMz-specCID.getPrecursorPeak().getMz()) >= 0.1f &&
						Math.abs(precursorMz-specETD.getPrecursorPeak().getMz()) >= 0.1f)
				{
					System.out.println(annotation);
					System.out.println(specCID.getTitle());
					System.out.println(specETD.getTitle());
					System.out.println(Math.abs(precursorMz-specCID.getPrecursorPeak().getMz()));
					System.out.println(Math.abs(precursorMz-specETD.getPrecursorPeak().getMz()));
				}
					
				assert(Math.abs(precursorMz-specCID.getPrecursorPeak().getMz()) < 0.1f ||
						Math.abs(precursorMz-specETD.getPrecursorPeak().getMz()) < 0.1f);
			}
			specCID.setAnnotation(new Peptide(annotation));
			specETD.setAnnotation(new Peptide(annotation));
			newContainerCID.add(specCID);
			newContainerETD.add(specETD);
		}
		newContainerCID.outputMgfFile(outputFileCID);
		newContainerETD.outputMgfFile(outputFileETD);
		System.out.println("Done");		
	}
	
	public static void makeAnnotatedSpectra() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/MSGF2D/rescored/TrypCID_0.txt";
		String specFileName1 = System.getProperty("user.home")+"/Research/MSGF2D/CID_paired_01_tryp.mgf";
		String specFileName2 = System.getProperty("user.home")+"/Research/MSGF2D/CID_paired_03_tryp.mgf";
		String outputFileName = System.getProperty("user.home")+"/Research/MSGF2D/annotatedSpectra/annotated_Tryp_CID.mgf";
		
		SpectraContainer newContainer = new SpectraContainer();
		float threshold = 2.8147016E-11f;	// TrypCID
//		float threshold = 1.1997423E-11f;	// TrypETD
//		float threshold = 1.0247307E-11f;	// LysNETD
//		float threshold = 1.8861999E-11f;	// LysNCID
		
//		cidThreshold = 1.8861999E-11f;	// LysNCID 3%
//		etdThreshold = 1.0247307E-11f;	// LysNETD 3%
//		cidThreshold = 2.8147016E-11f;	// TrypCID 3%
//		etdThreshold = 1.1997423E-11f;	// TrypETD 3%
		
		String s;
		Hashtable<String, String> pepTable = new Hashtable<String, String>();
		BufferedLineReader in = new BufferedLineReader(fileName);
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("SpecNum") || s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			if(token.length != 7)
				continue;
			float specProb = Float.parseFloat(token[6]);
			if(specProb < threshold)
			{
				String pep = token[4].substring(token[4].indexOf('.')+1, token[4].lastIndexOf('.'));
				if(pepTable.get(pep) == null)
					pepTable.put(pep, s);
				else
				{
					String[] token2 = pepTable.get(pep).split("\t");
					assert(token2.length == 7);
					float existingSpecProb = Float.parseFloat(token2[6]);
					if(specProb < existingSpecProb)
						pepTable.put(pep, s);
				}
			}
		}
		
		SpectraContainer container1 = new SpectraContainer(specFileName1, new MgfSpectrumParser());
		SpectraContainer container2 = new SpectraContainer(specFileName2, new MgfSpectrumParser());
		for(String str : pepTable.values())
		{
			String[] token = str.split("\t");
			int specNum = Integer.parseInt(token[0]);
			float precursorMz = Float.parseFloat(token[1]);
			String annotation = token[4].substring(token[4].indexOf('.')+1, token[4].lastIndexOf('.'));
			Spectrum spec1 = null;
			Spectrum spec2 = null;
			if(specNum-1 < container1.size())
				spec1 = container1.get(specNum-1);
			if(specNum-1 < container2.size())
				spec2 = container2.get(specNum-1);
			
			Spectrum spec = null;
			if(spec1 == null)
				spec = spec2;
			else if(spec2 == null)
				spec = spec1;
			else
			{
				if(Math.abs(precursorMz-spec1.getPrecursorPeak().getMz()) < Math.abs(precursorMz-spec2.getPrecursorPeak().getMz()))
					spec = spec1;
				else
					spec = spec2;
				assert(Math.abs(precursorMz-spec.getPrecursorPeak().getMz()) < 0.01f);
			}
			spec.setAnnotation(new Peptide(annotation));
			newContainer.add(spec);
		}
		newContainer.outputMgfFile(outputFileName);
		System.out.println("Done");
	}
	
	public static void vennDiagram() throws Exception
	{
		/*
		String fileNameCID = System.getProperty("user.home")+"/Research/MSGF2D/rescored/TrypCID_1.txt";
		String fileNameETD = System.getProperty("user.home")+"/Research/MSGF2D/rescored/TrypETD_1.txt";
		float cidThreshold = 4.85e-11f;
		float etdThreshold = 2.33e-11f;
//		float cidThreshold = 7.35e-12f;
//		float etdThreshold = 2.70e-12f;
		*/
		
		// 1.16136406E-10 11351, 9.211624E-11 6350
		// 6.785344E-11 6000, 4.3867715E-11 4524
		
		String enzyme = "Tryp";
		boolean isDecoy = false;
		String fileNameCID = System.getProperty("user.home")+"/Research/Data/HeckWhole/AnnotatedSpectra/annotatedHeckWholeSum_CID_"+enzyme+"CID_"+(isDecoy ? 1 : 0) + ".txt";
		String fileNameETD = System.getProperty("user.home")+"/Research/MSGF2D/rescored/"+enzyme+"ETD_"+(isDecoy ? 1 : 0) + ".txt";

		float cidThreshold, etdThreshold;
		if(enzyme.equalsIgnoreCase("LysN"))
		{
			cidThreshold = 1.8861999E-11f;	// LysNCID 3%
			etdThreshold = 1.0247307E-11f;	// LysNETD 3%
		}
		else
		{
			cidThreshold = 2.8147016E-11f;	// TrypCID 3%
			etdThreshold = 1.1997423E-11f;	// TrypETD 3%
		}
		
		//		float cidThreshold = 3.89e-11f;
//		float etdThreshold = 1.60e-11f;
//		float cidThreshold = 5.01e-12f;
//		float etdThreshold = 2.30e-12f;
		
		int numCID = 0;
		int numETD = 0;
		int numCIDOnly = 0;
		int numIntersection = 0;
		int numETDOnly = 0;
		int numConflict = 0;
		
		HashSet<String> cidOnly = new HashSet<String>();
		HashSet<String> etdOnly = new HashSet<String>();
		HashSet<String> both = new HashSet<String>();
		HashSet<String> confl = new HashSet<String>();
		
		Hashtable<Integer, ArrayList<String>> etdTable = new Hashtable<Integer, ArrayList<String>>();
		BufferedLineReader in = new BufferedLineReader(fileNameETD);
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("SpecNum") || s.startsWith("#"))
				continue;
			int specNum = Integer.parseInt(s.split("\t")[0]);
			ArrayList<String> list = etdTable.get(specNum);
			if(list == null)
			{
				ArrayList<String> newList = new ArrayList<String>();
				newList.add(s);
				etdTable.put(specNum, newList);
			}
			else
				list.add(s);
		}
		
		in = new BufferedLineReader(fileNameCID);
		
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("SpecNum") || s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			int specNum = Integer.parseInt(token[0]);
			ArrayList<String> etdResultList = etdTable.get(specNum);
			assert(etdResultList != null && etdResultList.size() >= 1 && etdResultList.size()<=2);
			String etdResult = null;
			if(etdResultList.size() == 1)
			{
				etdResult = etdResultList.get(0);
			}
			else if(etdResultList.size() == 2)
			{
				float mz0 = Float.parseFloat(etdResultList.get(0).split("\t")[1]);
				float mz1 = Float.parseFloat(etdResultList.get(1).split("\t")[1]);
				float mz = Float.parseFloat(token[1]);
				if(Math.abs(mz-mz0) < Math.abs(mz-mz1))
					etdResult = etdResultList.get(0);
				else
					etdResult = etdResultList.get(1);
			}
			String idCID = token[4];
			String idETD = etdResult.split("\t")[4];
			
			float specProbCID = Float.parseFloat(token[6]);
			float specProbETD = Float.parseFloat(etdResult.split("\t")[6]);
			
			if(specProbCID < cidThreshold)
				numCID++;
			if(specProbETD < etdThreshold)
				numETD++;
			if(specProbCID < cidThreshold && specProbETD < etdThreshold)
			{
				numIntersection++;
				if(!idCID.equalsIgnoreCase(idETD))
				{
					numConflict++;
					System.out.println(s);
					System.out.println(etdResult);
					confl.add(idCID+idETD);
				}
				else
					both.add(idCID);
			}
			else if(specProbCID < cidThreshold)
			{
				numCIDOnly++;
				cidOnly.add(idCID);
			}
			else if(specProbETD < etdThreshold)
			{
				numETDOnly++;
				etdOnly.add(idETD);
			}
		}			
		
		HashSet<String> cidExclusive = new HashSet<String>();
		HashSet<String> etdExclusive = new HashSet<String>();
		
		for(String id : cidOnly)
			if(!both.contains(id))
				cidExclusive.add(id);
		for(String id : etdOnly)
			if(!both.contains(id))
				etdExclusive.add(id);
		
		System.out.println("CID only: " + numCIDOnly);
		System.out.println("ETD only: " + numETDOnly);
		System.out.println("Intersection: " + numIntersection);
		System.out.println("Conflict: " + numConflict);
		System.out.println("NumCID: " + numCID);
		System.out.println("NumETD: " + numETD);
		
		System.out.println("CID only (peptide): " + cidExclusive.size());
		System.out.println("ETD only (peptide): " + etdExclusive.size());
		System.out.println("Intersection (peptide): " + both.size());
		System.out.println("Conflict (peptide): " + confl.size());
	}
	
	public static void best() throws Exception
	{
		/*
		String fileNameCID = System.getProperty("user.home")+"/Research/MSGF2D/rescored/LysNCID_0.txt";
		String fileNameETD = System.getProperty("user.home")+"/Research/MSGF2D/rescored/LysNETD_0.txt";
		String fileNameSum = System.getProperty("user.home")+"/Research/MSGF2D/rescored/LysNSum_0.txt";
		
		Hashtable<Integer, ArrayList<String>> etdTable = new Hashtable<Integer, ArrayList<String>>();
		BufferedLineReader in = new BufferedLineReader(fileNameETD);
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("SpecNum") || s.startsWith("#"))
				continue;
			int specNum = Integer.parseInt(s.split("\t")[0]);
			ArrayList<String> list = etdTable.get(specNum);
			if(list == null)
			{
				ArrayList<String> newList = new ArrayList<String>();
				newList.add(s);
				etdTable.put(specNum, newList);
			}
			else
				list.add(s);
		}
		
		in = new BufferedLineReader(fileNameCID);
		
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("SpecNum") || s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			int specNum = Integer.parseInt(token[0]);
			ArrayList<String> etdResultList = etdTable.get(specNum);
			assert(etdResultList != null && etdResultList.size() >= 1 && etdResultList.size()<=2);
			String etdResult = null;
			if(etdResultList.size() == 1)
			{
				etdResult = etdResultList.get(0);
			}
			else if(etdResultList.size() == 2)
			{
				float mz0 = Float.parseFloat(etdResultList.get(0).split("\t")[1]);
				float mz1 = Float.parseFloat(etdResultList.get(1).split("\t")[1]);
				float mz = Float.parseFloat(token[1]);
				if(Math.abs(mz-mz0) < Math.abs(mz-mz1))
					etdResult = etdResultList.get(0);
				else
					etdResult = etdResultList.get(1);
			}
			
			float specProbCID = Float.parseFloat(token[6]);
			float specProbETD = Float.parseFloat(etdResult.split("\t")[6]);
			if(specProbCID < specProbETD)
			else
		}		
		
		in = new BufferedLineReader(fileNameSum);
		
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("SpecNum") || s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			int specNum = Integer.parseInt(token[0]);
			ArrayList<String> etdResultList = etdTable.get(specNum);
			assert(etdResultList != null && etdResultList.size() >= 1 && etdResultList.size()<=2);
			String etdResult = null;
			if(etdResultList.size() == 1)
			{
				etdResult = etdResultList.get(0);
			}
			else if(etdResultList.size() == 2)
			{
				float mz0 = Float.parseFloat(etdResultList.get(0).split("\t")[1]);
				float mz1 = Float.parseFloat(etdResultList.get(1).split("\t")[1]);
				float mz = Float.parseFloat(token[1]);
				if(Math.abs(mz-mz0) < Math.abs(mz-mz1))
					etdResult = etdResultList.get(0);
				else
					etdResult = etdResultList.get(1);
			}
			
			float specProbCID = Float.parseFloat(token[6]);
			float specProbETD = Float.parseFloat(etdResult.split("\t")[6]);
			if(specProbCID < specProbETD)
				System.out.println(s);
			else
				System.out.println(etdResult);
		}		
		*/
		
	}	
	
	public static void better() throws Exception
	{
		String fileNameCID = System.getProperty("user.home")+"/Research/MSGF2D/rescored/LysNCID_0.txt";
		String fileNameETD = System.getProperty("user.home")+"/Research/MSGF2D/rescored/LysNETD_0.txt";
		
		Hashtable<Integer, ArrayList<String>> etdTable = new Hashtable<Integer, ArrayList<String>>();
		BufferedLineReader in = new BufferedLineReader(fileNameETD);
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("SpecNum") || s.startsWith("#"))
				continue;
			int specNum = Integer.parseInt(s.split("\t")[0]);
			ArrayList<String> list = etdTable.get(specNum);
			if(list == null)
			{
				ArrayList<String> newList = new ArrayList<String>();
				newList.add(s);
				etdTable.put(specNum, newList);
			}
			else
				list.add(s);
		}
		
		in = new BufferedLineReader(fileNameCID);
		
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("SpecNum") || s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			int specNum = Integer.parseInt(token[0]);
			ArrayList<String> etdResultList = etdTable.get(specNum);
			assert(etdResultList != null && etdResultList.size() >= 1 && etdResultList.size()<=2);
			String etdResult = null;
			if(etdResultList.size() == 1)
			{
				etdResult = etdResultList.get(0);
			}
			else if(etdResultList.size() == 2)
			{
				float mz0 = Float.parseFloat(etdResultList.get(0).split("\t")[1]);
				float mz1 = Float.parseFloat(etdResultList.get(1).split("\t")[1]);
				float mz = Float.parseFloat(token[1]);
				if(Math.abs(mz-mz0) < Math.abs(mz-mz1))
					etdResult = etdResultList.get(0);
				else
					etdResult = etdResultList.get(1);
			}
			
			float specProbCID = Float.parseFloat(token[6]);
			float specProbETD = Float.parseFloat(etdResult.split("\t")[6]);
			if(specProbCID < specProbETD)
				System.out.println(s);
			else
				System.out.println(etdResult);
		}		
	}
	
	public static void intersection() throws Exception
	{
		int numPrecursorMismatch = 0;
		String fileNameCID = System.getProperty("user.home")+"/Research/MSGF2D/rescored/TrypCID_0.txt";
		String fileNameETD = System.getProperty("user.home")+"/Research/MSGF2D/rescored/TrypETD_0.txt";
		
		Hashtable<String, String> etdTable = new Hashtable<String, String>();
		BufferedLineReader in = new BufferedLineReader(fileNameETD);
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("SpecNum") || s.startsWith("#"))
				continue;
			int specNum = Integer.parseInt(s.split("\t")[0]);
			etdTable.put(specNum+":"+s.split("\t")[1], s);
		}
		
		in = new BufferedLineReader(fileNameCID);
		
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("SpecNum") || s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			int specNum = Integer.parseInt(token[0]);
			String etdResult = etdTable.get(specNum+":"+token[1]);
			if(etdResult == null)
			{
				numPrecursorMismatch++;
				continue;
			}
			String idCID = s.split("\t")[4];
			String idETD = etdResult.split("\t")[4];
			if(idCID.equalsIgnoreCase(idETD))
				System.out.println(s+"\t"+etdResult.split("\t")[6]);
		}
//		System.out.println("NumPrecursorMismatch: " + numPrecursorMismatch);
	}
	

	public static void mergeGRIDResults() throws Exception
	{
		String dirName = System.getProperty("user.home")+"/Research/MSGF2D/paired";
		File dir = new File(dirName);
		if(!dir.isDirectory())
			return;
		ArrayList<Integer> scanNumList = new ArrayList<Integer>();
		Hashtable<Integer, String> result = new Hashtable<Integer, String>();
		
		boolean labelPrinted = false;
		for(File f : dir.listFiles())
		{
			if(!f.getName().endsWith("1_0.txt"))
				continue;
			BufferedLineReader in = new BufferedLineReader(f.getPath());
			String s = in.readLine();
			if(labelPrinted == false)
			{
				System.out.println(s);
				labelPrinted = true;
			}
			while((s=in.readLine()) != null)
			{
				String[] token = s.split("\t");
				int scanNum = Integer.parseInt(token[1]);
				scanNumList.add(scanNum);
				result.put(scanNum, s);
			}
		}
		Collections.sort(scanNumList);
		for(Integer scanNum : scanNumList)
			System.out.println(result.get(scanNum));
	}
}
