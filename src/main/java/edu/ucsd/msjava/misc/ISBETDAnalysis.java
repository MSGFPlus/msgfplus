package edu.ucsd.msjava.misc;

import java.io.File;
import java.util.ArrayList;

import edu.ucsd.msjava.parser.PSM;
import edu.ucsd.msjava.parser.PSMList;
import edu.ucsd.msjava.parser.PepXMLParser;


public class ISBETDAnalysis {
	public static void main(String argv[]) throws Exception
	{
		generateAnnotatedSpectra();
//		generateISBResults();
	}
	
	public static void generateISBResults() throws Exception
	{
//		String dirName = System.getProperty("user.home")+"/Research/Data/ISBETD/BioRep2TechRep1_ETD/iprophet/";
//		String dirName = System.getProperty("user.home")+"/Research/Data/ISBETD/BioRep2TechRep1_ETD/SEQ_YeastCombNR_20070207_ForwDecoy/";
		String dirName = System.getProperty("user.home")+"/Research/Data/ISBETD/BioRep2TechRep1_ETD/OMScp_YeastCombNR_20070207_ForwDecoy/";
		String fileName = "interact-ipro.pep.xml";
		String pepXMLFileName = dirName+"/"+fileName;
		PSMList<PSM> psmList = PepXMLParser.parse(pepXMLFileName);
		ArrayList<String> scoreNames = null;
		if(psmList.size() > 0)
		{
			scoreNames = psmList.get(0).getScoreNames();
		}
		System.out.print("#SpecFile\tScanNum\tCharge\tAnnotation");
		for(String name : scoreNames)
			System.out.print("\t"+name);
		System.out.println("\tProtein");
		for(PSM psm : psmList)
		{
//			if(psm.getPeptide().isModified())
//				continue;
			if(!psm.getProtein().contains("DECOY"))
			{
				System.out.print(psm.getSpecFileName()+"\t"+psm.getScanNum()+"\t"+psm.getCharge()+"\t"+
						psm.getPrecedingResidue()+"."+psm.getPeptideStr()+"."+psm.getSucceedingResidue());
				for(String name : scoreNames)
				{
					System.out.print("\t"+psm.getScore(name));
						
				}
				System.out.println("\t"+psm.getProtein());
			}
		}
	}
	
	public static void generateAnnotatedSpectra()
	{
		String pepXMLFileName = System.getProperty("user.home")+"/Research/Data/ISBETD/BioRep2TechRep1_ETD/iprophet/interact-ipro.pep.xml";
		PSMList<PSM> psmList = PepXMLParser.parse(pepXMLFileName).getDistinctivePeptideSet();
		String outputFileName = System.getProperty("user.home")+"/Research/Data/ISBETD/AnnotatedSpectra/annotatedISBETD_ipro09.mgf";
		psmList.outputMgf(new File(System.getProperty("user.home")+"/Research/Data/ISBETD/BioRep2TechRep1_ETD"), new File(outputFileName), "interprophet", 0.9f, true);
	}
	
}
