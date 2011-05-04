package misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import msutil.SpectraIterator;
import msutil.Spectrum;

import parser.BufferedLineReader;
import parser.MgfSpectrumParser;
import parser.MzXMLSpectraIterator;
import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;

public class IPRGStudy {
	public static void main(String argv[]) throws Exception
	{
//		areAllSpectraProcessed();
//		extractModResults();
//		count();
//		enzymeCleavageTest();
//		convertIntoIPRGReport();
//		generateIPRGReport();
		mergeHeckAndiPRGSpectra();
	}

	public static void mergeHeckAndiPRGSpectra() throws Exception
	{
		String heckSpecName = System.getProperty("user.home")+"/Research/Data/HeckRevision/AnnotatedSpectra/ETD_Tryp_Confident.mgf";
		String iprgSpecName = System.getProperty("user.home")+"/Research/Data/ABRF/2011/allFrxns/annotatedABRFETD.mgf";
		String outputFileName = System.getProperty("user.home")+"/Research/Data/ABRF/2011/allFrxns/annotatedABRFETDWithHeckCharge2.mgf";
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFileName)));
		
		SpectraIterator itr = new SpectraIterator(heckSpecName, new MgfSpectrumParser());
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			String pepStr = spec.getAnnotationStr();
			if(pepStr.endsWith("K"))
				spec.outputMgf(out);
		}
		itr = new SpectraIterator(iprgSpecName, new MgfSpectrumParser());
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			spec.outputMgf(out);
		}
		out.flush();
		out.close();
		System.out.println("Done");
	}
	
	public static void convertIntoIPRGReport() throws Exception
	{
		HashMap<Character,String> ptmMap = new HashMap<Character,String>();
		ptmMap.put('m', "Oxidation");
		ptmMap.put('k', "Acetylation");
//		ptmMap.put('e', "Pyro-glu");
//		ptmMap.put('q', "Pyro-glu");
//		ptmMap.put('n', "Deamidation");
//		ptmMap.put('q', "Deamidation");
		convertIntoIPRGReport(
				"/home/sangtaekim/Research/Data/ABRF/StudyFiles/MSGFDB_1209_30ppm_OxAcetyl.txt", ptmMap,
				"/home/sangtaekim/Research/Data/ABRF/StudyFiles/FinalReport_OxAcetyl.txt", false);
	}
	
	public static void convertIntoIPRGReport(String msgfOutput, HashMap<Character,String> ptmMap, String iPRGOutput, boolean nTerm) throws Exception
	{
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(iPRGOutput)));
		
		String s;
		BufferedLineReader in = new BufferedLineReader(msgfOutput);
		in.readLine();	// header
		out.println("Scan number\tPrecursor m/z\tMass error (ppm)\tPrecursor charge\tPeptide Sequence\tModifications\tProtein Accession(s)\tSpectral Probability\tBetter than 1% FDR threshod?");
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < 11)
				continue;
			int scanNum = Integer.parseInt(token[1]);
			float precursorMz = Float.parseFloat(token[3]);
			int charge = Integer.parseInt(token[5]);
			float precursorError = Float.parseFloat(token[4])/charge;
			String annotation = token[6];
			String pepStr = annotation.substring(annotation.indexOf('.')+1, annotation.lastIndexOf('.'));
			StringBuffer modifications = new StringBuffer();
			
			boolean isPTMOK = true;
			for(int i=0; i<pepStr.length(); i++)
			{
				char residue = pepStr.charAt(i);
				if(Character.isLowerCase(residue))
				{
					if(nTerm && residue != 'm' && i != 0)
					{
						isPTMOK = false;
						break;
					}
					String ptm = ptmMap.get(residue);
					modifications.append((i+1)+","+Character.toUpperCase(residue)+","+ptm+";");
				}
			}
			
			if(!isPTMOK)
				continue;
			String prot = token[7];
			float specProb = Float.parseFloat(token[10]);
			out.println(scanNum+"\t"+precursorMz+"\t"+precursorError+"\t"+charge+"\t"+pepStr+"\t"+modifications.toString()+"\t"+prot+"\t"+specProb);
		}
		in.close();
		out.close();
		System.out.println("Done");
	}
	
	public static void generateIPRGReport() throws Exception
	{
		String mergedResults = System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/FinalReport_Merged.txt";

		File temp1 = File.createTempFile("temp", "merge");
		temp1.deleteOnExit();
		mergeSearchResults(
				System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/FinalReport_Nomod.txt",
				System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/FinalReport_OxDeamd.txt",
				temp1.getPath());

		File temp2 = File.createTempFile("temp", "merge");
		temp2.deleteOnExit();
		mergeSearchResults(
				temp1.getPath(),
				System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/FinalReport_OxPyroglu.txt",
				temp2.getPath());
		
		mergeSearchResults(
				temp2.getPath(),
				System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/FinalReport_OxAcetyl.txt",
				mergedResults);
		
	}
	
	public static void mergeSearchResults(String targetResults, String decoyResults, String mergedResults) throws Exception
	{
//		String targetResults = System.getProperty("user.home")+"/Research/Data/ISBControl/Mix_7/ORBITRAP/MSGFDB/MSGFDB_S6.txt";
//		String decoyResults = System.getProperty("user.home")+"/Research/Data/ISBControl/Mix_7/ORBITRAP/MSGFDB/MSGFDB_S10.txt";
//		String mergedResults = System.getProperty("user.home")+"/Research/Data/ISBControl/Mix_7/ORBITRAP/MSGFDB/MSGFDB_S24.txt";
		int specProbColumn = 7;
		int keyColumn = 0;
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
			if(targetResult == null)
			{
				out.println(s);
				continue;
			}
			else
			{
				float targetSpecProb = Float.parseFloat(targetResult.split("\t")[specProbColumn]);
				float decoySpecProb = Float.parseFloat(token[specProbColumn]);
				if(targetSpecProb <= decoySpecProb)
					out.println(targetResult);
				else if(targetSpecProb > decoySpecProb)
					out.println(s);
			}
		}		
		
		in.close();
		out.close();
		System.out.println("Done");
	}	
	public static void enzymeCleavageTest() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/MSGFDB_1209_30ppm_Merged_FDR01.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		in.readLine();	// header
		String s;
		int[] numSpecs = new int[100];
		int[] numSpecsPreceedingK = new int[100];
		int[] numSpecsEndingK = new int[100];
		
		int numAllSpecs = 0;
		int numAllSpecsPreceedingK = 0;
		int numAllSpecsEndingK = 0;
		
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			String annotationStr = token[6];
			int charge = Integer.parseInt(token[5]);
			numSpecs[charge]++;
			numAllSpecs++;
			if(annotationStr.charAt(0) == 'K')
			{
				numSpecsPreceedingK[charge]++;
				numAllSpecsPreceedingK++;
			}
			if(annotationStr.charAt(annotationStr.lastIndexOf('.')-1) == 'K')
			{
				numSpecsEndingK[charge]++;
				numAllSpecsEndingK++;
			}
		}
		for(int c=2; c<=10; c++)
			System.out.println(c+"\t"+numSpecs[c]+"\t"+numSpecsPreceedingK[c]/(float)numSpecs[c]+"\t"+numSpecsEndingK[c]/(float)numSpecs[c]);
		System.out.println("All\t"+numAllSpecs+"\t"+numAllSpecsPreceedingK/(float)numAllSpecs+"\t"+numAllSpecsEndingK/(float)numAllSpecs);
	}
	public static void count() throws Exception
	{
		String fileName1 = System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/ModOnly.txt";
		String fileName2 = System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/NomodOnly.txt";
		HashSet<String> scanNumSet = new HashSet<String>();
		BufferedLineReader in = new BufferedLineReader(fileName1);
		in.readLine();	// header
		String s;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			scanNumSet.add(token[1]);
		}
		in = new BufferedLineReader(fileName2);
		in.readLine();	// header
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			scanNumSet.add(token[1]);
		}
		System.out.println(scanNumSet.size());
	}
	
	public static void extractModResults() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/ABRF/StudyFiles/MSGFDB_1207_30ppm_Mod.txt";
		BufferedLineReader in = new BufferedLineReader(fileName);
		System.out.println(in.readLine());	// header
		String s;
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length < 6)
				continue;
			String pepStr = token[6].substring(token[6].indexOf('.')+1, token[6].lastIndexOf('.'));
			if(!pepStr.toUpperCase().equals(pepStr))
				System.out.println(s);
		}
		in.close();
	}
	
	public static void areAllSpectraProcessed() throws Exception
	{
		String fileName = "/home/sangtaekim/Research/Data/ABRF/StudyFiles/MSGFDB_1101_PMErr_20ppm_Deamd.txt";
		fileName = "/home/sangtaekim/Developments/MS_Java/bin/msgfdb_highPMtol_reversed.raw.1.txt";
		String specFileName = "/home/sangtaekim/Research/Data/ABRF/StudyFiles/D100930_yeast_SCX10S_rak_ft8E_pc_01.mzXML";
		specFileName = "/home/sangtaekim/Research/ToolDistribution/MSGFDBTest/ABRF_Ab_HC_DTT_NIPIA_ArgC_37C_On_102909.mzXML";
		
		String s;
		BufferedLineReader in = new BufferedLineReader(fileName);
		in.readLine();
		HashSet<Integer> scanNumSet = new HashSet<Integer>();
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			int scanNum = Integer.parseInt(token[1]);
			scanNumSet.add(scanNum);
		}
		System.out.println("NumProcessedSpecs\t"+scanNumSet.size());
		
		MzXMLSpectraIterator itr = new MzXMLSpectraIterator(specFileName);
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			int scanNum = spec.getScanNum();
			if(!scanNumSet.contains(scanNum))
				System.out.println("Missing\t"+scanNum);
		}
	}
}

