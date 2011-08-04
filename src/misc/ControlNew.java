package misc;

import java.io.*;
import java.util.Arrays;

import parser.MzXMLSpectraMap;
import parser.PSM;
import parser.PSMList;
import parser.PepXMLParser;

import msutil.SpectraContainer;
import msutil.Spectrum;

public class ControlNew {
	public static void main(String argv[])
	{
		if(argv.length != 1)
			printUsageAndExit();
		File dir = new File(argv[0]);
		if(!dir.isDirectory())
			printUsageAndExit();
		else
			generateAnnotatedSpectra(dir);

		System.out.println("Done"); 
	}
	
	public static void printUsageAndExit()
	{
		System.out.println("usage: java ControlNew directory");
		System.exit(0);
	}
	
	
	public static String[] controlProteins = {
		"ATBOG",
		"P00634",
		"AMY_BACLI",
		"P00711",
		"P02666",
		"P00722",
		"P02754",
		"P00921",
		"P00432",
		"P00006",
		"P46406",
		"P00489",
		"P00946",
		"P02188",
		"P02602",
		"P01012",
		"Q29443",
		"P02769"
	};
	
	public static void generateAnnotatedSpectra(File root)
	{
		File mzXMLDir = null;
		File searchResultsDir = null;
		for(File f : root.listFiles())
		{
			if(f.getName().equalsIgnoreCase("mzXML"))
				mzXMLDir = f;
			else if(f.getName().equalsIgnoreCase("Search_Results"))
				searchResultsDir = f;
		}
		
		if(mzXMLDir == null || searchResultsDir == null || !mzXMLDir.isDirectory() || !searchResultsDir.isDirectory())
		{
			System.out.println("Wrong directory!");
			System.exit(0);
		}
		
		File[] mzXMLFiles = mzXMLDir.listFiles(new FileFilter.FileExtFilter("mzXML"));
		File[] pepProphetFiles = searchResultsDir.listFiles(new FileFilter.FilePrefixFilter("interact"));
		
		if(mzXMLFiles.length != pepProphetFiles.length)
		{
			System.out.println("Number of spectra and annotation files are different");
			System.exit(0);
		}

		Arrays.sort(mzXMLFiles);
		Arrays.sort(pepProphetFiles);

		String annoDirName = "Annotated_Specs";
		File annotatedSpecDir = new File(root, annoDirName);
		annotatedSpecDir.mkdir();

		for(int i=0; i<mzXMLFiles.length; i++)
		{
			String fileName = mzXMLFiles[i].getName();
			String name = fileName.substring(0, fileName.lastIndexOf('.'));
			System.out.println("Processing " + name);
			String resultFileName = root.getPath() + File.separatorChar + annoDirName + File.separatorChar + name + "_Anno.mgf";
			getAnnotatedSpectra(mzXMLFiles[i], pepProphetFiles[i]).outputMgfFile(resultFileName);
		}
	}
	
	public static SpectraContainer getAnnotatedSpectra(File mzXML, File pepProphetXML)
	{
		MzXMLSpectraMap specDB = new MzXMLSpectraMap(mzXML.getPath());
		SpectraContainer specListWithAnnotation = new SpectraContainer();

		PSMList<PSM> psmList = PepXMLParser.parse(pepProphetXML.getPath());
		for(PSM psm : psmList)
		{
			if(psm.getPeptideStr() == null)
				continue;
			boolean controlProtein = false;
			for(String keyword : controlProteins)
			{
				if(psm.getProtein().contains(keyword))
				{
					controlProtein = true;
					break;
				}
			}
			if(!controlProtein)
				continue;
			Spectrum spec = specDB.getSpectrumBySpecIndex(psm.getScanNum());
			spec.setAnnotation(psm.getPeptide());
			spec.setCharge(psm.getCharge());
			spec.setScanNum(psm.getScanNum());
			spec.setTitle("Prob_" + psm.getScore("peptideprophet"));
			specListWithAnnotation.add(spec);
		}
		return specListWithAnnotation;
	}
}
