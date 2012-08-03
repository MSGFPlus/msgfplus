package edu.ucsd.msjava.misc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;

import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.SpectraContainer;
import edu.ucsd.msjava.msutil.SpectraMap;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.SpectrumAccessorBySpecIndex;
import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.parser.MgfSpectrumParser;
import edu.ucsd.msjava.parser.MzXMLSpectraMap;

public class AnnotatedSpecGenerator {
	public static void main(String argv[]) throws Exception
	{
		File resultFile = null;
		File specDir = null;
		File outputFile = null;
		int specIndexCol = -1;
		int specFileCol = -1;
		int peptideCol = -1;
		int chargeCol = -1;
		int scoreCol = -1;
		boolean isGreaterBetter = false;
		float threshold = Float.MAX_VALUE;
		boolean uniquePeptide = false;
		boolean hasHeader = true;
//		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		int i=0;
		while(i<argv.length)
		{
     		if(argv[i].equalsIgnoreCase("-r"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				resultFile = new File(argv[i+1]);
				if(!resultFile.exists())
					printUsageAndExit(argv[i+1] + " doesn't exist.");
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-d"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				specDir = new File(argv[i+1]);
				if(!specDir.exists() || !specDir.isDirectory())
					printUsageAndExit(argv[i+1] + " doesn't exist or not a directory!");
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-o"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				outputFile = new File(argv[i+1]);
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-n"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				specIndexCol = Integer.parseInt(argv[i+1]);
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-f"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				specFileCol = Integer.parseInt(argv[i+1]);
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-p"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				peptideCol = Integer.parseInt(argv[i+1]);
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-c"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				chargeCol = Integer.parseInt(argv[i+1]);
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-s"))
			{
				if(i+2 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				scoreCol = Integer.parseInt(argv[i+1]);
				if(argv[i+2].equalsIgnoreCase("1"))
					isGreaterBetter = true;
				i += 3;
			}
     		else if(argv[i].equalsIgnoreCase("-t"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				threshold = Float.parseFloat(argv[i+1]);
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-u"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				if(argv[i+1].equalsIgnoreCase("1"))
					uniquePeptide = true;
				i += 2;
			}
     		else if(argv[i].equalsIgnoreCase("-h"))
			{
				if(i+1 >= argv.length)
					printUsageAndExit("Illegal parameter: " + argv[i]);
				if(argv[i+1].equalsIgnoreCase("0"))
					hasHeader = false;
				i += 2;
			}
//			else if(argv[i].equalsIgnoreCase("-fixMod"))
//			{
//				// 0: No mod, 1: Carbamidomethyl C, 2: Carboxymethyl C
//				if(argv[i+1].equalsIgnoreCase("0"))
//					aaSet = AminoAcidSet.getStandardAminoAcidSet();
//				else if(argv[i+1].equalsIgnoreCase("1"))
//					aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
//				else if(argv[i+1].equalsIgnoreCase("2"))
//					aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarboxymethylatedCys();
//				else
//					printUsageAndExit("Illigal -fixMod parameter: " + argv[i+1]);
//				i += 2;
//			}     		
     		else
     			printUsageAndExit("Illegal parameter!");
		}
		
		if(resultFile == null)
			printUsageAndExit("-r resultFileName is missing");
		if(outputFile == null)
			printUsageAndExit("-o outputFileName is missing");
		if(specDir == null)
			printUsageAndExit("-d specDir is missing");
		if(specIndexCol < 0)
			printUsageAndExit("-n scanNumCol is missing");
		if(specFileCol < 0)
			printUsageAndExit("-f specFileCol is missing");
		if(peptideCol < 0)
			printUsageAndExit("-p peptideCol is missing");
		if(chargeCol < 0)
			printUsageAndExit("-c chargeCol is missing");
		if(scoreCol < 0)
			printUsageAndExit("-s scoreCol 0/1 is missing");
		
		if(threshold == Float.MAX_VALUE)
		{
			if(isGreaterBetter)
				threshold = Float.MIN_VALUE;
			else
				threshold = Float.MAX_VALUE;
		}
		
		generateAnnotatedSpectra(resultFile, specDir, outputFile, specIndexCol, specFileCol, peptideCol, chargeCol, scoreCol, isGreaterBetter,
				threshold, uniquePeptide, hasHeader);
		
	}
	
	public static void printUsageAndExit(String message)
	{
		System.out.println(message);
		System.out.println("usage: java AnnotatedSpecGenerator\n" +
				"\t-r resultFileName\n" +
				"\t-d specDir\n" +
				"\t-o outputFileName\n" +
				"\t-f specFileCol\n" +
				"\t-n specIndexCol\n" +
				"\t-p peptideCol\n" +
				"\t-c chargeCol\n" +
				"\t-s scoreCol 0/1 (0: smaller is better, 1: larger is better)\n" +
				"\t[-t threshold] \n" +
				"\t[-u 0/1] (0: one spectrum per peptide, 1: no restriction (default))\n" +
				"\t[-h 0/1] (0: no header, 1: header (default))\n"
//				"\t[-fixMod 0/1/2] (0: NoCysteineProtection, 1: CarbamidomethyC (default), 2: CarboxymethylC)\n"
				);
		System.exit(-1);
	}
	
	public static void generateAnnotatedSpectra(File resultFile, File specDir, File outputFile, int specIndexCol, int specFileCol, int pepCol, int chargeCol,
			int scoreCol, boolean isGreaterBetter, float threshold, boolean uniquePeptide, boolean hasHeader) throws Exception
	{
		String s;
		Hashtable<String, String> pepTable = null;
		ArrayList<String> resultList = null;
		if(uniquePeptide)
			pepTable = new Hashtable<String, String>();
		else
			resultList = new ArrayList<String>();
		
		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		
		BufferedLineReader in = new BufferedLineReader(resultFile.getPath());
		if(hasHeader)
			in.readLine();
		
		while((s=in.readLine()) != null)
		{
			String[] token = s.split("\t");
			if(token.length <= specIndexCol || token.length <= specFileCol || token.length <= pepCol || token.length <= scoreCol)
				continue;
			String pep = token[pepCol];
//			if(pep.matches("[A-Z]\\..+\\.[A-Z]"))
			if(pep.matches(".\\..+\\.."))
				pep = pep.substring(pep.indexOf('.')+1, pep.lastIndexOf('.'));
			float score = Float.parseFloat(token[scoreCol]);
			
			if(!isGreaterBetter && score < threshold || isGreaterBetter && score > threshold)
			{
				if(uniquePeptide)
				{
					if(pepTable.get(pep) == null)
						pepTable.put(pep, s);
					else
					{
						String[] token2 = pepTable.get(pep).split("\t");
						float existingScore = Float.parseFloat(token2[scoreCol]);
						if(score < existingScore)
							pepTable.put(pep, s);
					}
				}
				else
					resultList.add(s);
			}
		}
		
		Iterator<String> itr;
		if(uniquePeptide)
			itr = pepTable.values().iterator();
		else
			itr = resultList.iterator();

		HashMap<String,SpectrumAccessorBySpecIndex> specAccessorMap = new HashMap<String,SpectrumAccessorBySpecIndex>(); 
		while(itr.hasNext())
		{
			String str = itr.next();
			String[] token = str.split("\t");
			
			String pep = token[pepCol];
//			if(pep.matches(".*\\.[A-Z]+\\..*"))
			if(pep.matches(".\\..+\\.."))
				pep = pep.substring(pep.indexOf('.')+1, pep.lastIndexOf('.'));
//			if(!pep.contains("79.966"))
//				continue;
//			if(!pep.contains("phos"))
//				continue;
			
			String specFileName = token[specFileCol];
			specFileName = new File(specFileName).getName();
			SpectrumAccessorBySpecIndex specMap = specAccessorMap.get(specFileName);
			if(specMap == null)
			{
				String ext = specFileName.substring(specFileName.lastIndexOf('.'));
				if(ext.equalsIgnoreCase(".mzXML"))
					specMap = new MzXMLSpectraMap(specDir.getPath()+File.separator+specFileName);
				else if(ext.equalsIgnoreCase(".mgf"))
					specMap = new SpectraMap(specDir.getPath()+File.separator+specFileName, new MgfSpectrumParser());
				else
				{
					System.out.println("Unrecognized spectrum format: " + specFileName);
					System.exit(-1);
				}
				specAccessorMap.put(specFileName, specMap);
			}
			
			int specIndex = Integer.parseInt(token[specIndexCol]);
			int charge = Integer.parseInt(token[chargeCol]);
			Spectrum spec = specMap.getSpectrumBySpecIndex(specIndex);
			if(spec == null)
			{
				System.out.println(specFileName+":"+specIndex+" is not available!");
				System.exit(-1);
			}
			
			out.println("BEGIN IONS");
			out.println("TITLE=" + specFileName+":"+specIndex+" MSLevel="+spec.getMSLevel());
			out.println("SEQ=" + pep);
			float precursorMz = spec.getPrecursorPeak().getMz();
			out.println("PEPMASS=" + precursorMz);
			if(spec.getScanNum() > 0)
			out.println("SCANS=" + spec.getScanNum());
			out.println("CHARGE=" + charge+ (charge>0 ? "+" : ""));
			for(Peak p : spec)
				if(p.getIntensity() > 0)
					out.println(p.getMz() + "\t" + p.getIntensity());
			out.println("END IONS");
		}
		out.close();
		System.out.println("Done");		
	}
	
}
