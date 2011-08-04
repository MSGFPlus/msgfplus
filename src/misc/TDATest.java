package misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import msgf.Tolerance;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.SpectraIterator;
import msutil.SpectraMap;
import msutil.Spectrum;

import fdr.Pair;
import fdr.TargetDecoyPSMSet;

import parser.BufferedLineReader;
import parser.MgfSpectrumParser;
import parser.SortedSpectraIterator;

import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;



public class TDATest {
	public static void main(String argv[]) throws Exception
	{
//		checkDBSize();
//			summaryOMSSA(i, true);
//		make100ShuffledDB();
//		calcDBSize();
//		reverseFilteredDB();
//		for(int i=0; i<=99; i++)
//			summarizeSearch19(i, false);
//		filterTwoPassResults();
//		checkDBSize();
//		dbSizeComparison(7);
//		mergeSearchResults();
//		System.out.println(generateDetailedResults(21,false));
//		computeEfficiency();
//		summarize();
//		splitByCharges();
//		mergeSearchResults();
		postProcessProtein();
	}

	public static void postProcessProtein() throws Exception
	{
		String dbFileName = "/home/sangtaekim/Research/Data/ISBControl/Mix_7/ORBITRAP/Database/18mix.fasta";
		SuffixArray sa = new SuffixArray(new SuffixArraySequence(dbFileName));
		
		String resultName = "/home/sangtaekim/Research/Data/TDATest/ASMS/TDATest_S3N_Target.txt";
		int pepColumn = 6;
		int protColumn = 7;

		BufferedLineReader in = new BufferedLineReader(resultName);
		
		System.out.println(in.readLine());
		String s;
		while((s=in.readLine()) != null)
		{
			String token[] = s.split("\t");
			if(token.length <= protColumn)
				continue;
			String annotation = token[pepColumn];
			String protein = token[protColumn];
			String str = annotation.replaceAll("[\\._]", "");
			
			if(protein.startsWith("\"Y") || protein.startsWith("Y"))
			{
				for(String protMatch : sa.getAllMatchingAnnotations(str))
				{
					if(!protMatch.startsWith("Y"))
						protein = protMatch.split("\\s+")[0];
				}
			}
			
			System.out.print(token[0]);
			for(int i=1; i<token.length; i++)
			{
				if(i != protColumn)
					System.out.print("\t"+token[i]);
				else
					System.out.print("\t"+protein);
			}
			System.out.println();
		}
	}
	
	static int parentMassErrorThresholdPPM = 30;
	private static float fdrThreshold = 0.01f;
	static float[] efficiency;
	
	static {
		efficiency = new float[40];
		float effMassOnly = 0.291932565f;
		float effMassProt = 0.002190685f;
		float effSpecMass = 0.1308691f;
		float effAll = 0.001238213f;
		
		for(int i=0; i<efficiency.length; i++)
			efficiency[i] = effAll;
		
		efficiency[1] = effMassOnly;
		efficiency[3] = effMassOnly;
		efficiency[28] = effMassOnly;
		efficiency[29] = effMassOnly;
		
		efficiency[2] = effMassProt; 
		efficiency[4] = effMassProt; 
		
		efficiency[5] = effSpecMass;
		efficiency[30] = efficiency[32] = 1-0.7157701f;
		efficiency[31] = efficiency[33] = 1-0.93468535f;
	}
	
	public static void summarize() throws Exception
	{
		boolean isPeptideLevel = false;
		for(int i=34; i<=35; i++)
		{
			if(i != 19 && i != 10 && i != 14 && i != 26)
			{
				String output = summaryMSGF(i, isPeptideLevel);
				if(output != null)
					System.out.println(output);
				if(i==21)
				{
					String output6 = summaryMSGF(6, isPeptideLevel);
					String[] token6 = output6.split("\t");
					String[] token = output.split("\t");
					System.out.print("6+21");
					for(int j=1; j<=9; j++)
						System.out.print("\t"+(Integer.parseInt(token6[j])+Integer.parseInt(token[j])));
					float factualFDR;
					int num;
					if(!isPeptideLevel)
						num = Integer.parseInt(token6[1])+Integer.parseInt(token[1]);
					else
						num = Integer.parseInt(token6[3])+Integer.parseInt(token[3]);
					factualFDR = (num-Integer.parseInt(token6[8])-Integer.parseInt(token[8]))/(float)num;
					float adjFactFDR = factualFDR*Float.parseFloat(token[11])/Float.parseFloat(token[10]);
					System.out.println("\t"+factualFDR+"\t"+adjFactFDR);
				}
			}
		}
	}
	
	public static void splitByCharges() throws Exception
	{
		for(int i=1; i<=33; i++)
		{
			if(i != 19 && i != 10 && i != 12 && i != 14 && i != 26 && i != 22 && i != 23 && i != 26)
			{
				System.out.println(i);
				summaryMSGFByCharges(i);
				System.out.println();
			}
		}
	}	
	public static void computeEfficiency() throws Exception
	{
		int expNum = 31;
		String dirName = System.getProperty("user.home")+"/Research/Data/TDATest";
		String fileName = dirName+"/MSGFDB_S"+expNum+"_Decoy.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		in.readLine();
		String s;
		int numPSMs = 0;
		int spec = 0;
		int prot = 0;
		int mass = 0;
		int specProt = 0;
		int specMass = 0;
		int protMass = 0;
		int all = 0;
		
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < 13)
				continue;
			numPSMs++;
			
			float minError = Float.MAX_VALUE;
			float precursorMz = Float.parseFloat(token[4]);
			float pmerror = Float.parseFloat(token[5]);
			int charge = Integer.parseInt(token[6]);
			for(int err=-1; err<=2; err++)
			{
				float parentMass = (precursorMz-err*(float)(Composition.C13-Composition.C)-(float)Composition.PROTON)*charge;
				float error = (pmerror-err*(float)Composition.PROTON)/parentMass*1e6f;
				if(Math.abs(error) < Math.abs(minError))
					minError = error;
			}	
			float adjError = minError;
			
			int specFilter=0, protFilter=0, massFilter=0, allFilters=0;
			if(token[2].startsWith("ABRF"))
				specFilter = 1;
			if(token[8].startsWith("REV_Y") || token[8].startsWith("REV_AT"))
				protFilter = 1;
			if(Math.abs(adjError) >= parentMassErrorThresholdPPM)
				massFilter = 1;
			if(specFilter+protFilter+massFilter > 0)
				allFilters = 1;
			spec += specFilter;
			prot += protFilter;
			mass += massFilter;
			specProt += 1-Math.max(1-specFilter-protFilter, 0);
			protMass += 1-Math.max(1-protFilter-massFilter, 0);
			specMass += 1-Math.max(1-specFilter-massFilter, 0);
			all += allFilters;
		}
		System.out.println("Spec\t"+(spec/(float)numPSMs));
		System.out.println("Prot\t"+(prot/(float)numPSMs));
		System.out.println("Mass\t"+(mass/(float)numPSMs));
		System.out.println("SpecProt\t"+(specProt/(float)numPSMs));
		System.out.println("ProtMass\t"+(protMass/(float)numPSMs));
		System.out.println("SpecMass\t"+(specMass/(float)numPSMs));
		System.out.println("All\t"+(all/(float)numPSMs));
	}
	
	public static int generateDetailedResults(int expNum, boolean printResult) throws Exception
	{
		String dirName = System.getProperty("user.home")+"/Research/Data/TDATest";
		String fileName = dirName+"/MSGFDB_S"+expNum+".txt";
		
		StringBuffer output = new StringBuffer();
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		s = in.readLine();	// header
		String[] token = s.split("\t");
		output.append("Rank");
		for(int i=0; i<token.length; i++)
		{
			output.append("\t"+token[i]);
			if(i==5)
				output.append("\tAdjErr");
			else if(i==8)
				output.append("\tIsDecoy\tNumDecoy\tNumTarget\tFDR");
		}
		output.append("\tSpecFilter\tProtFilter\tMassFilter\tAllFilters\tFactFDR\tAdjFactFDR");
		output.append("\n");
		boolean isErrorPPM = false;
		if(s.split("\t")[5].contains("ppm"))
			isErrorPPM = true;
		int rank = 0;
		int numDecoy = 0;
		int numTarget = 0;
		int numFactualPSMs = 0;
		int numNonFactualPSMs = 0;
		
		while((s=in.readLine()) != null)
		{
			token = s.split("\t");
			if(token.length < 13)
				continue;
			rank++;
			int isDecoy = 0;
			String prot = token[8];
			if(prot.startsWith("REV") || prot.startsWith("SHFL") || prot.startsWith("XXX"))
				isDecoy = 1;
			
			if(isDecoy == 0)
				output.append(rank);
			float adjError = Float.parseFloat(token[5]);
			for(int i=0; i<token.length; i++)
			{
				if(isDecoy == 0)
					output.append("\t"+token[i]);
				if(i==5)	// PMError
				{
					if(!isErrorPPM)
					{
						float minError = Float.MAX_VALUE;
						float precursorMz = Float.parseFloat(token[4]);
						float pmerror = Float.parseFloat(token[5]);
						int charge = Integer.parseInt(token[6]);
						for(int err=-1; err<=2; err++)
						{
							float parentMass = (precursorMz-err*(float)(Composition.C13-Composition.C)-(float)Composition.PROTON)*charge;
							float error = (pmerror-err*(float)Composition.PROTON)/parentMass*1e6f;
							if(Math.abs(error) < Math.abs(minError))
								minError = error;
						}	
						adjError = minError;
					}
					if(isDecoy == 0)
						output.append("\t"+adjError);
				}
				else if(i==8)
				{
					numDecoy += isDecoy;
					numTarget = rank-numDecoy;
					float fdr = numDecoy/(float)numTarget;
					if(isDecoy == 0)
						output.append("\t"+isDecoy+"\t"+numDecoy+"\t"+numTarget+"\t"+fdr);
				}
			}
			if(isDecoy == 0)
			{
				int specFilter=0, protFilter=0, massFilter=0, allFilters=0;
				if(token[2].startsWith("ABRF"))
					specFilter = 1;
				if(token[8].startsWith("Y") || token[8].startsWith("AT"))
					protFilter = 1;
				if(Math.abs(adjError) >= parentMassErrorThresholdPPM)
					massFilter = 1;
				if(specFilter+protFilter+massFilter > 0)
					allFilters = 1;
				if(allFilters == 1)
					numNonFactualPSMs++;
				else
					numFactualPSMs++;
				output.append("\t"+specFilter+"\t"+protFilter+"\t"+massFilter+"\t"+allFilters+"\t"+(numNonFactualPSMs/(float)numTarget)
						+"\t"+(numNonFactualPSMs/(float)numTarget)/(1-efficiency[expNum]));
				output.append("\n");
			}
		}
		String results = output.toString();
		String[] lines = results.split("\n");
		output = new StringBuffer();
		float minFDR = Float.MAX_VALUE;
		float minFactFDR = Float.MAX_VALUE;
		for(int lineNum = lines.length-1; lineNum > 0; lineNum--)
		{
			token = lines[lineNum].split("\t");
			for(int i=0; i<token.length; i++)
			{
				if(i==14)	// FDR
				{
					float fdr = Float.parseFloat(token[14]);
					if(fdr > minFDR)
						fdr = minFDR;
					else
						minFDR = fdr;
					output.append(fdr+"\t");
				}
				else if(i==24)	// FactFDR
				{
					float factFDR = Float.parseFloat(token[24]);
					if(factFDR > minFactFDR)
						factFDR = minFactFDR;
					else
						minFactFDR = factFDR;
					output.append(factFDR+"\t"+factFDR/(1-efficiency[expNum])+"\n");
				}
				else if(i!=25)
					output.append(token[i]+"\t");
			}
		}
		output.append(lines[0]+"\n");
		lines = output.toString().split("\n");
		float adjFactFDR = 0;
		int numPSM = 0;
		for(int lineNum=lines.length-1; lineNum>=0; lineNum--)
		{
			if(lineNum < lines.length-1)
			{
				token = lines[lineNum].split("\t");
				adjFactFDR = Float.parseFloat(token[25]);
				numPSM = Integer.parseInt(token[13]);
				if(!printResult && adjFactFDR > fdrThreshold)
					break;
			}
			if(printResult)
				System.out.println(lines[lineNum]);
		}
		return numPSM-1;
	}
	
	public static void dbSizeComparison(int expNum) throws Exception
	{
		System.out.println("Best\tSecondBest\tDecoyBest\tBest-Second\tBest-Decoy\tRatio");
		String dirName = System.getProperty("user.home")+"/Research/Data/ISBControl/Mix_7/ORBITRAP/MSGFDB0105";
		String targetFileName = dirName+"/MSGFDB_S"+expNum+"_Target.txt";
		String decoyFileName = dirName+"/MSGFDB_S"+expNum+"_Decoy.txt";
		String s;
		HashMap<Integer, ArrayList<String>> target = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> decoy = new HashMap<Integer, ArrayList<String>>();
		
		// read target
		BufferedLineReader in = new BufferedLineReader(targetFileName);
		in.readLine();	// header
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < 11)
				continue;
			int scanNum = Integer.parseInt(token[1]);
			ArrayList<String> resList = target.get(scanNum);
			if(resList == null)
				resList = new ArrayList<String>();
			resList.add(s);
			target.put(scanNum, resList);
		}
		
		// read decoy
		in = new BufferedLineReader(decoyFileName);
		in.readLine();	// header
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < 11)
				continue;
			int scanNum = Integer.parseInt(token[1]);
			ArrayList<String> resList = decoy.get(scanNum);
			if(resList == null)
				resList = new ArrayList<String>();
			resList.add(s);
			decoy.put(scanNum, resList);
		}		
		
		float threshold = 1;
		if(expNum == 5)
			threshold = 3.521849E-7f;
		else if(expNum == 7)
			threshold = 2.5315636E-10f;
		else if(expNum == 16)
			threshold = 2.9182195E-11f;
		
		for(int scanNum : target.keySet())
		{
			ArrayList<String> list = target.get(scanNum);
			assert(list != null): "Null List: " + scanNum;
			if(list.size() < 2 || decoy.get(scanNum) == null)
			{
				System.out.println("1: " + scanNum);
				continue;
			}
			String best = list.get(0);
			String[] token = best.split("\t");
			float bestPSMSpecProb = Float.parseFloat(token[11]);
			if(bestPSMSpecProb > threshold)
			{
				continue;
			}
			
			float bestMSGFScore = -(float)Math.log10(Float.parseFloat(token[11]));
			float secondBestMSGFScore = bestMSGFScore;
			for(int i=1; i<list.size(); i++)
			{
				float score = -(float)Math.log10(Float.parseFloat(list.get(i).split("\t")[11]));
				if(score < bestMSGFScore)
				{
					secondBestMSGFScore = score;
					break;
				}
			}
			if(secondBestMSGFScore >= bestMSGFScore)
			{
				System.out.println("3: " + scanNum + "\t" + bestMSGFScore + "\t" + secondBestMSGFScore);
				assert(false);
				continue;
			}
			String decoyBest = decoy.get(scanNum).get(0);
			// Positive PSM
			float decoyBestMSGFScore = -(float)Math.log10(Float.parseFloat(decoyBest.split("\t")[11]));
			if(decoyBestMSGFScore >= bestMSGFScore)
				continue;
			System.out.format("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
					bestMSGFScore,secondBestMSGFScore,decoyBestMSGFScore,
					bestMSGFScore-secondBestMSGFScore,
					bestMSGFScore-decoyBestMSGFScore,
					(bestMSGFScore-secondBestMSGFScore)/(float)(bestMSGFScore-decoyBestMSGFScore));
		}
	}
	
	public static void calcNumDistinctPeptides() throws Exception
	{
		int length = 7;
		String dbFileName = "/home/sangtaekim/Research/Data/ISBControl/Mix_7/ORBITRAP/Database/18mix.fasta";
		SuffixArray sa = new SuffixArray(new SuffixArraySequence(dbFileName));
		
	}
	
	public static void filterTwoPassResults() throws Exception
	{
		String resultPass1 = System.getProperty("user.home")+"/Research/Data/TDATest/MSGFDB_S6_FDR01.txt";
		String resultPass2 = System.getProperty("user.home")+"/Research/Data/TDATest/MSGFDB_S20.txt";
		String filteredPass1 = System.getProperty("user.home")+"/Research/Data/TDATest/MSGFDB_S21.txt";

//		String resultPass1 = System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/MSGFDB_1207_30ppm_Nomod_FDR001.txt";;
//		String resultPass2 = System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/MSGFDB_1207_30ppm_Mod.txt";;
//		String filteredPass1 = System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/MSGFDB_1207_30ppm_TwoPass.txt";;
		int keyColumn = 2;
		
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(filteredPass1)));
		
		HashMap<String,String> targetMap = new HashMap<String,String>();
		String s;
		BufferedLineReader in = new BufferedLineReader(resultPass1);
		in.readLine();	// header
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			String title = token[keyColumn];
			targetMap.put(title, s);
		}
		in.close();
		
		in = new BufferedLineReader(resultPass2);
		out.println(in.readLine());
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			String title = token[keyColumn];
			String targetResult = targetMap.get(title);
			if(targetResult == null)
			{
				out.println(s);
			}
//			else
//				out.println(targetResult);
		}		
		
		in.close();
		out.close();
		System.out.println("Done");
	}
	
	public static void mergeSearchResults() throws Exception
	{
//		String mergedResults = System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/MSGFDB_1209_30ppm_Merged.txt";
//
//		File temp1 = File.createTempFile("temp", "merge");
//		temp1.deleteOnExit();
//		mergeSearchResults(
//				System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/MSGFDB_1209_30ppm_Nomod.txt",
//				System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/MSGFDB_1209_30ppm_OxDeamd.txt",
//				temp1.getPath());
//
//		File temp2 = File.createTempFile("temp", "merge");
//		temp2.deleteOnExit();
//		mergeSearchResults(
//				temp1.getPath(),
//				System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/MSGFDB_1209_30ppm_OxPyroglu.txt",
//				temp2.getPath());
//		
//		mergeSearchResults(
//				temp2.getPath(),
//				System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/MSGFDB_1209_30ppm_OxAcetyl.txt",
//				mergedResults);

		String targetResults = System.getProperty("user.home")+"/Research/Data/TDATest/MSGFDB_S34_Target.txt";
		String decoyResults = System.getProperty("user.home")+"/Research/Data/TDATest/MSGFDB_S34_Decoy.txt";
		String mergedResults = System.getProperty("user.home")+"/Research/Data/TDATest/MSGFDB_S35.txt";
		mergeSearchResults(targetResults, decoyResults, mergedResults);
	}
	
	public static void mergeSearchResults(String targetResults, String decoyResults, String mergedResults) throws Exception
	{
//		int specProbColumn = 10;
//		int keyColumn = 1;
		int specProbColumn = 11;
		int keyColumn = 2;
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(mergedResults)));
		
		HashMap<String,String> targetMap = new HashMap<String,String>();
		String s;
		BufferedLineReader in = new BufferedLineReader(targetResults);
		in.readLine();	// header
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < keyColumn || token.length < specProbColumn)
				continue;
			String title = token[keyColumn];
			targetMap.put(title, s);
		}
		in.close();
		
		in = new BufferedLineReader(decoyResults);
		out.println(in.readLine());
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < keyColumn || token.length < specProbColumn)
				continue;
			String title = token[keyColumn];
			String targetResult = targetMap.get(title);
			StringBuffer decoyResult = new StringBuffer();
			for(int i=0; i<token.length-1; i++)
			{
//				if(i != 8)
					decoyResult.append(token[i]+"\t");
//				else
//					decoyResult.append("REV_"+token[i]+"\t");
			}
			decoyResult.append(token[token.length-1]);
			if(targetResult == null)
			{
				out.println(decoyResult);
				continue;
			}
			else
			{
				float targetSpecProb = Float.parseFloat(targetResult.split("\t")[specProbColumn]);
				float decoySpecProb = Float.parseFloat(token[specProbColumn]);
				if(targetSpecProb <= decoySpecProb)
					out.println(targetResult);
				else if(targetSpecProb > decoySpecProb)
					out.println(decoyResult);
			}
		}		
		
		in.close();
		out.close();
		System.out.println("Done");
	}
	
	public static void summarizeSearch19(int copyNum, boolean isPeptideLevel) throws Exception
	{
		File dir = new File(System.getProperty("user.home")+"/Research/Data/TDATest/MSGFDB_S19");
		String decoyPrefix = "SHFL";
		
		TargetDecoyPSMSet psmSet;
		
		psmSet = new TargetDecoyPSMSet(
			new File(dir.getPath()+File.separator+"MSGFDB_S19_"+copyNum+".txt"), 
			"\t", 
			true,
			11, 
			false, 
			0,
			1, 
			7,
			null,
			8, 
			decoyPrefix
			);


		File tempFile = File.createTempFile("FDRTest", "Temp");
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(tempFile)));
		float fdrThreshold = 0.01f;
		if(!isPeptideLevel)
			psmSet.writeResults(out, fdrThreshold, 1f, true);
		else
			psmSet.writeResults(out, 1, fdrThreshold, true);
		ArrayList<Number> numPSMs = processTempFileMSGF(tempFile, fdrThreshold, isPeptideLevel);
		
		float factualFDR = (numPSMs.get(0).intValue()-numPSMs.get(numPSMs.size()-1).intValue())/(float)numPSMs.get(0).intValue();
		System.out.print(copyNum+"\t"+numPSMs.get(0));
		for(int i=1; i<numPSMs.size(); i++)
			System.out.print("\t"+numPSMs.get(i));
		System.out.println("\t"+factualFDR);
	}
	
	public static void reverseFilteredDB() throws Exception
	{
		String fileName = "/home/sangtaekim/Research/Data/ISBControl/Mix_7/ORBITRAP/Database/ISBYeast_Rev_Comb_OMSPass2.fasta";
		String outFileName = "/home/sangtaekim/Research/Data/ISBControl/Mix_7/ORBITRAP/Database/ISBYeast_Rev_Comb_OMSPass2_RevComb.fasta";
		BufferedLineReader in = new BufferedLineReader(fileName);
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outFileName)));
		String s;
		boolean reverse = false;
		while((s = in.readLine()) != null)
		{
			if(s.startsWith(">"))
			{
				if(s.startsWith(">REV"))
					reverse = true;
				else
				{
					reverse = false;
					out.println(s);
				}
			}
			else if(reverse == false)
				out.println(s);
		}
		
		in = new BufferedLineReader(fileName);
		
		StringBuffer protein = null;
		String annotation = null;
		while((s = in.readLine()) != null)
		{
			if(s.startsWith(">"))	// start of a protein
			{
				if(annotation != null)
				{
					StringBuffer rev = new StringBuffer();
					for(int i=protein.length()-1; i>=0; i--)
						rev.append(protein.charAt(i));
					out.println(">REV_" + annotation);
					out.println(rev.toString().trim());
				}
				if(s.startsWith(">REV"))
				{
					reverse = true;
					annotation = null;
				}
				else
				{
					reverse = false;
					annotation = s.substring(1);
					protein = new StringBuffer();
				}
			}
			else if(!reverse)
				protein.append(s);
		}
		if(protein != null && annotation != null)
		{
			StringBuffer rev = new StringBuffer();
			for(int i=protein.length()-1; i>=0; i--)
				rev.append(protein.charAt(i));
			out.println(">REV_" + annotation);
			out.println(rev.toString());
		}
		in.close();
		out.close();
		System.out.println("Done");
	}
	
	public static void calcDBSize() throws Exception
	{
		int length = 20;
		String fileName = "/home/sangtaekim/Research/Data/ISBControl/Mix_7/ORBITRAP/Database/ISB_Rev_Comb.fasta";
//		fileName = "/home/sangtaekim/Research/Data/SProt/uniprot_sprot.fasta";
		BufferedLineReader in = new BufferedLineReader(fileName);
		String s;
		int numPep = 0;
		int protSize = -1;
		int numAA = 0;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith(">"))
			{
				if(protSize > 0)
					numPep += protSize-length+1;
				protSize = 0;
				continue;
			}
			else
			{
				protSize += s.length();
				numAA += s.length();
			}
		}
		if(protSize > 0)
			numPep += protSize-length+1;
		
		System.out.println("NumPeptides\t"+numPep);
		System.out.println("NumAA\t"+numAA);
		System.out.println("Ratio\t"+numPep/(float)numAA);
	}
	
	public static void make100ShuffledDB() throws Exception
	{
		File dbFile = new File("/home/sangtaekim/Research/Data/ISBControl/Mix_7/ORBITRAP/Database/18mix.fasta");
		for(int i=0; i<100; i++)
		{
			File outputFile = new File("/home/sangtaekim/Research/Data/ISBControl/Mix_7/ORBITRAP/Database/18mixShuffles/18mix_shuffle_"+i+".fasta");
			msdbsearch.ShuffleDB.shuffleDB(dbFile.getPath(), outputFile.getPath(), true);
		}
	}
	
	public static void summaryOMSSA(int expNum, boolean isPeptideLevel) throws Exception
	{
		File dir = new File(System.getProperty("user.home")+"/Research/OMSSA/omssa-2.1.9.linux");
		String decoyPrefix;
		
		if(expNum == 3 || expNum == 4)
			decoyPrefix = "SHFL";
		else
			decoyPrefix = "REV";
		
		TargetDecoyPSMSet psmSet;
		
		if(expNum != 7)
		{
			psmSet = new TargetDecoyPSMSet(
				new File(dir.getPath()+File.separator+"OMSSA_S"+expNum+".csv"), 
				",", 
				true,
				3, 
				false, 
				1,
				0, 
				2,
				null,
				9, 
				decoyPrefix
				);
		}
		else
		{
			psmSet = new TargetDecoyPSMSet(
					new File(dir.getPath()+File.separator+"OMSSA_S"+expNum+"_Target.csv"), 
					new File(dir.getPath()+File.separator+"OMSSA_S"+expNum+"_Decoy.csv"), 
					",", 
					true,
					3, 
					false, 
					1,
					0, 
					2,
					null
				);
		}

		File tempFile = File.createTempFile("FDRTest", "Temp");
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(tempFile)));
		float fdrThreshold = 0.01f;
		if(!isPeptideLevel)
			psmSet.writeResults(out, fdrThreshold, 1f, true);
		else
			psmSet.writeResults(out, 1, fdrThreshold, true);
		ArrayList<Number> numPSMs = processTempFileOMSSA(tempFile, fdrThreshold, isPeptideLevel);
		
		tempFile.deleteOnExit();
		System.out.print(numPSMs.get(0));
		for(int i=1; i<numPSMs.size(); i++)
			System.out.print("\t"+numPSMs.get(i));
		System.out.println();
	}	
	
	public static String summaryMSGF(int expNum, boolean isPeptideLevel) throws Exception
	{
		return summaryMSGF(expNum, isPeptideLevel, null);
	}
	
	// minChrage,maxCharge are inclusive
	public static void summaryMSGFByCharges(int expNum) throws Exception
	{
		boolean isPeptideLevel = false;
//		System.out.println("Experiment\t"+expNum);
		
		StringBuffer output = new StringBuffer();
		int chargeColumn = 6;
		int[][] numbers = new int[4][9];
		for(int charge=2; charge<=4; charge++)
		{
			ArrayList<Pair<Integer,ArrayList<String>>> reqStrList = new ArrayList<Pair<Integer,ArrayList<String>>>();;
			ArrayList<String> list = new ArrayList<String>();
			if(charge == 2 || charge == 3)
				list.add(String.valueOf(charge));
			else
			{
				for(int c=4; c<=100; c++)
					list.add(String.valueOf(c));
			}
			reqStrList.add(new Pair<Integer,ArrayList<String>>(chargeColumn, list));
			String result = summaryMSGF(expNum, isPeptideLevel, reqStrList);
			String[] token = result.split("\t");
			for(int i=1; i<10; i++)
			{
				int num = Integer.parseInt(token[i]);
				numbers[charge-2][i-1] += num;
			}
		}
		
		for(int i=0; i<3; i++)
			for(int j=0; j<9; j++)
				numbers[3][j] += numbers[i][j];
		
		for(int charge=2; charge<=5; charge++)
		{
			if(charge==2 || charge == 3)
				System.out.print("C"+charge);
			else if(charge == 4)
				System.out.print("C4+");
			else
				System.out.print("Sum");
			for(int i=0; i<9; i++)
			{
				if(i==2 || i==3)
					continue;
				System.out.print("\t"+numbers[charge-2][i]);
			}
			float factualFDR;
			if(!isPeptideLevel)
				factualFDR = (numbers[charge-2][0]-numbers[charge-2][7])/(float)numbers[charge-2][0];
			else
				factualFDR = (numbers[charge-2][2]-numbers[charge-2][7])/(float)numbers[charge-2][2];
			System.out.print("\t"+factualFDR);
			System.out.print("\t"+factualFDR/(1-efficiency[expNum]));
			System.out.println();
		}
	}
	
	public static String summaryMSGF(int expNum, boolean isPeptideLevel, ArrayList<Pair<Integer,ArrayList<String>>> reqStrList) throws Exception
	{
		File dir = new File(System.getProperty("user.home")+"/Research/Data/TDATest/");
//		File dir = new File(System.getProperty("user.home")+"/Research/Data/ISBControl/Mix_7/ORBITRAP/MSGFDB");
		String decoyPrefix;
		
		if(expNum == 3 || expNum == 4 || expNum == 26 || expNum == 27 || expNum == 29 || expNum == 34 || expNum == 35)
			decoyPrefix = "SHFL";
		else
			decoyPrefix = "REV";
		
		TargetDecoyPSMSet psmSet;
		
		if(expNum != 7 && expNum != 16 && expNum != 18 && expNum != 25 && expNum != 26 && expNum != 28 && expNum != 29 && expNum != 30 && expNum != 31
				&& expNum != 34)
		{
			File file = new File(dir.getPath()+File.separator+"MSGFDB_S"+expNum+".txt"); 
			if(!file.exists())
				return null;
			psmSet = new TargetDecoyPSMSet(
				file, 
				"\t", 
				true,
				11, 
				false, 
				0,
				1, 
				7,
				reqStrList,
				8, 
				decoyPrefix
				);
		}
		else
		{
			File targetFile = new File(dir.getPath()+File.separator+"MSGFDB_S"+expNum+"_Target.txt");
			File decoyFile = new File(dir.getPath()+File.separator+"MSGFDB_S"+expNum+"_Decoy.txt");
			if(!targetFile.exists() || !decoyFile.exists())
				return null;
			psmSet = new TargetDecoyPSMSet(
					targetFile, 
					decoyFile, 
					"\t", 
					true,
					11, 
					false,
					0,
					1, 
					7,
					reqStrList,
					1
				);
		}

		float thresholdScore = psmSet.getThresholdScore(fdrThreshold, isPeptideLevel);
		File tempFile = File.createTempFile("FDRTest", "Temp");
		tempFile.deleteOnExit();
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(tempFile)));
		if(!isPeptideLevel)
			psmSet.writeResults(out, 1f, 1f, thresholdScore, true);
		else
			psmSet.writeResults(out, 1f, 1f, thresholdScore, true);
		
		// #targetPSMs #decoyPSMs #targetPeptides #decoyPeptides MassFilter ProtFilter SpecFilter AllFilters Unmod
		StringBuffer output = new StringBuffer();
		ArrayList<Number> numbers = processTempFileMSGF(tempFile, fdrThreshold, isPeptideLevel);
//		System.out.print(expNum);
		output.append(expNum);
		for(int i=0; i<numbers.size(); i++)
			output.append("\t"+numbers.get(i));
//			System.out.print("\t"+numbers.get(i));
		float factualFDR;
		if(!isPeptideLevel)
			factualFDR = (numbers.get(0).intValue()-numbers.get(numbers.size()-1).intValue())/(float)numbers.get(0).intValue();
		else
			factualFDR = (numbers.get(2).intValue()-numbers.get(numbers.size()-1).intValue())/(float)numbers.get(2).intValue();

		output.append("\t"+factualFDR);
//		System.out.print("\t"+factualFDR);
		output.append("\t"+factualFDR/(1-efficiency[expNum]));
//		System.out.print("\t"+factualFDR/(1-weight));
//		System.out.println();
		return output.toString();
	}
	
	static 	String[] control18 = {
		"P62739",
		"P00634",
		"P06278",
		"P00711",
		"P02666",
		"P00722",
		"P02754",
		"P00921",
		"P00432",
		"P62894",
		"P46406",
		"P00489",
		"P00946",
		"P68082",
		"P02602",
		"P01012",
		"Q29443",
		"P02769"
	};

	public static ArrayList<Number> processTempFileMSGF(File tempFile, float fdrThreshold, boolean usePeptideLevel) throws Exception
	{
		ArrayList<Number> numberList = new ArrayList<Number>();
		
		// FDR threhold: 0.01
		int numPSMs = 0;
		int numDecoyPSMs = 0;
		int numPSMWithErr = 0;
		int numPSMWithCorrProt = 0;
		int numPSMWithCorrSpec = 0;
		int numPassAll = 0;
		int numPSMsNoMod = 0;
		
		HashSet<String> pepSet = new HashSet<String>();
		HashSet<String> decoyPepSet = new HashSet<String>();
		HashSet<String> pspSetWithErr = new HashSet<String>();
		HashSet<String> pspSetWithCorrProt = new HashSet<String>();
		HashSet<String> pspSetWithCorrSpec = new HashSet<String>();
		HashSet<String> pspSetWithPassAll = new HashSet<String>();
		HashSet<String> pspSetWithNoMod = new HashSet<String>();
		
		BufferedLineReader in = new BufferedLineReader(tempFile.getPath());
		
		String header = in.readLine();	// header
		boolean isErrorDa;
		if(header.split("\t")[5].contains("Da"))
			isErrorDa = true;
		else
			isErrorDa = false;
		
		String s;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < 12)
				continue;
			String specFileName = token[2].substring(0, token[2].lastIndexOf(':'));
			int scanNum = Integer.parseInt(token[2].substring(token[2].lastIndexOf(':')+1));
			float precursorMz = Float.parseFloat(token[4]);
			float pmerror = Float.parseFloat(token[5]);
			int charge = Integer.parseInt(token[6]);
			String annotation = token[7];
			String pepStr = annotation.substring(annotation.indexOf('.')+1, annotation.lastIndexOf('.'));
			String prot = token[8];
//			float specProb = Float.parseFloat(token[11]);
//			float fdr = Float.parseFloat(token[14]);
//			float pepFDR = Float.parseFloat(token[15]);

			if(prot.startsWith("REV") || prot.startsWith("SHFL") || prot.startsWith("IPI:REV") || prot.startsWith("XXX"))
			{
				decoyPepSet.add(pepStr);
				numDecoyPSMs++;
				continue;
			}
			
			numPSMs++;
			pepSet.add(pepStr.toUpperCase());

			boolean passErr = false;
			// error check
			if(!isErrorDa)
			{
				pmerror = pmerror*((precursorMz-(float)Composition.ISOTOPE)*charge)/1e6f;
			}
			for(int err=-1; err<=2; err++)
			{
				float parentMass = (precursorMz-err*(float)(Composition.C13-Composition.C)-(float)Composition.PROTON)*charge;
				float adjErr = (pmerror-err*(float)Composition.PROTON)/parentMass*1e6f;
				if(Math.abs(adjErr) < parentMassErrorThresholdPPM)
				{
					numPSMWithErr++;
					pspSetWithErr.add(pepStr.toUpperCase());
					passErr = true;
					break;
				}
			}

			boolean passProt = false;
			// protein check
			boolean correctProt = false;
			for(String corr : control18)
			{
				if(prot.contains(corr))
				{
					correctProt = true;
					break;
				}
			}
			if(correctProt)
			{
				numPSMWithCorrProt++;
				pspSetWithCorrProt.add(pepStr.toUpperCase());
				passProt = true;
			}
			
			boolean passSpec = false;
			// spectrum check
			if(specFileName.startsWith("OR20070924"))
			{
				numPSMWithCorrSpec++;
				pspSetWithCorrSpec.add(pepStr.toUpperCase());
				passSpec = true;
			}
			
			if(passErr && passProt && passSpec)
			{
				numPassAll++;
				pspSetWithPassAll.add(pepStr.toUpperCase());

				// count unmodified peptides
				if(pepStr.toUpperCase().equals(pepStr))
				{
					numPSMsNoMod++;
					pspSetWithNoMod.add(pepStr.toUpperCase());
				}
			}
			else
			{
//				System.out.println("Debug");
			}
			
		}
		
		numberList.add(numPSMs);
		numberList.add(numDecoyPSMs);
		numberList.add(pepSet.size());	// number of target peptides
		numberList.add(decoyPepSet.size());	// number of decoy peptides
		
		if(!usePeptideLevel)
		{
			numberList.add(numPSMWithErr);
			numberList.add(numPSMWithCorrProt);
			numberList.add(numPSMWithCorrSpec);
			numberList.add(numPassAll);
			numberList.add(numPSMsNoMod);
		}
		else
		{
			numberList.add(pspSetWithErr.size());
			numberList.add(pspSetWithCorrProt.size());
			numberList.add(pspSetWithCorrSpec.size());
			numberList.add(pspSetWithPassAll.size());
			numberList.add(pspSetWithNoMod.size());
		}
		
		return numberList;
	}
	
	public static ArrayList<Number> processTempFileOMSSA(File tempFile, float fdrThreshold, boolean usePeptideLevel) throws Exception
	{
		ArrayList<Number> numberList = new ArrayList<Number>();
		
		// FDR threhold: 0.01
		int numPSMs = 0;
		int numDecoyPSMs = 0;
		int numPSMWithErr = 0;
		int numPSMWithCorrProt = 0;
		int numPSMWithCorrSpec = 0;
		int numPassAll = 0;
		int numPSMsNoMod = 0;
		
		HashSet<String> pepSet = new HashSet<String>();
		HashSet<String> decoyPepSet = new HashSet<String>();
		HashSet<String> pspSetWithErr = new HashSet<String>();
		HashSet<String> pspSetWithCorrProt = new HashSet<String>();
		HashSet<String> pspSetWithCorrSpec = new HashSet<String>();
		HashSet<String> pspSetWithPassAll = new HashSet<String>();
		HashSet<String> pspSetWithNoMod = new HashSet<String>();
		
		BufferedLineReader in = new BufferedLineReader(tempFile.getPath());
		
		String header = in.readLine();	// header
		
		String s;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split(",");
			if(token.length < 17)
				continue;
			
			String specFileName = token[1].substring(0, token[1].lastIndexOf(':'));
			int scanNum = Integer.parseInt(token[1].substring(token[1].lastIndexOf(':')+1));
			float spectrumPepMass = Float.parseFloat(token[4]);
			float theoPepMass = Float.parseFloat(token[token.length-5]);
			
			String pepStr = token[2];
			String prot = token[9];
			float eValue = Float.parseFloat(token[3]);
			float fdr = Float.parseFloat(token[token.length-2]);
			float pepFDR = Float.parseFloat(token[token.length-1]);

			if(prot.startsWith("REV") || prot.startsWith("SHFL"))
			{
				decoyPepSet.add(pepStr);
				numDecoyPSMs++;
				continue;
			}
			
			numPSMs++;
			pepSet.add(pepStr.toUpperCase());

			boolean passErr = false;
			// error check
			for(int err=-2; err<=2; err++)
			{
				float adjErr = ((spectrumPepMass-err*(float)Composition.PROTON)-theoPepMass)/theoPepMass*1e6f;
				if(Math.abs(adjErr) < parentMassErrorThresholdPPM)
				{
					numPSMWithErr++;
					pspSetWithErr.add(pepStr.toUpperCase());
					passErr = true;
					break;
				}
			}

			boolean passProt = false;
			// protein check
			boolean correctProt = false;
			for(String corr : control18)
			{
				if(prot.contains(corr))
				{
					correctProt = true;
					break;
				}
			}
			if(correctProt)
			{
				numPSMWithCorrProt++;
				pspSetWithCorrProt.add(pepStr.toUpperCase());
				passProt = true;
			}
			
			boolean passSpec = false;
			// spectrum check
			if(specFileName.startsWith("OR20070924"))
			{
				numPSMWithCorrSpec++;
				pspSetWithCorrSpec.add(pepStr.toUpperCase());
				passSpec = true;
			}
			
			if(passErr && passProt && passSpec)
			{
				numPassAll++;
				pspSetWithPassAll.add(pepStr.toUpperCase());

				// count unmodified peptides
				if(pepStr.toUpperCase().equals(pepStr))
				{
					numPSMsNoMod++;
					pspSetWithNoMod.add(pepStr.toUpperCase());
				}
			}
		}
		
		numberList.add(numPSMs);
		numberList.add(numDecoyPSMs);
		numberList.add(pepSet.size());	// number of target peptides
		numberList.add(decoyPepSet.size());	// number of decoy peptides
		
		if(!usePeptideLevel)
		{
			numberList.add(numPSMWithErr);
			numberList.add(numPSMWithCorrProt);
			numberList.add(numPSMWithCorrSpec);
			numberList.add(numPassAll);
			numberList.add(numPSMsNoMod);
		}
		else
		{
			numberList.add(pspSetWithErr.size());
			numberList.add(pspSetWithCorrProt.size());
			numberList.add(pspSetWithCorrSpec.size());
			numberList.add(pspSetWithPassAll.size());
			numberList.add(pspSetWithNoMod.size());
		}
		
		return numberList;
	}	
}
