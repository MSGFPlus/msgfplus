package edu.ucsd.msjava.misc;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;

import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.parser.InsPecTPSM;
import edu.ucsd.msjava.parser.InsPecTParser;
import edu.ucsd.msjava.parser.PSMList;


public class CompareSearchResults {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 2)
			printUsageAndExit();
		compare(new File(argv[0]), new File(argv[1]));
	}
	
	public static void printUsageAndExit()
	{
		System.out.println("Usage: java CompareSearchResults resultFile1 resultFile2");
		System.exit(-1);
	}
	
	public static void compare(File result1, File result2) throws Exception
	{
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();

		InsPecTParser parser1 = new InsPecTParser(aaSet);
		parser1.parse(result1.getPath());
		PSMList<InsPecTPSM> psmList1 = parser1.getPSMList("FDR", 0.01f, false);
		
		InsPecTParser parser2 = new InsPecTParser(aaSet);
		parser2.parse(result2.getPath());
		PSMList<InsPecTPSM> psmList2 = parser2.getPSMList("FDR", 0.01f, false);
		
		HashMap<Integer, InsPecTPSM> specIndexMap1 = new HashMap<Integer, InsPecTPSM>();
		for(InsPecTPSM psm : psmList1)
			specIndexMap1.put(psm.getScanNum()+1, psm);
		
		HashMap<Integer, InsPecTPSM> specIndexMap2 = new HashMap<Integer, InsPecTPSM>();
		for(InsPecTPSM psm : psmList2)
			specIndexMap2.put(psm.getSpecIndex(), psm);
		
		HashSet<Integer> common = new HashSet<Integer>();
		for(Integer specIndex : specIndexMap1.keySet())
		{
			if(specIndexMap2.containsKey(specIndex))
				common.add(specIndex);
		}
		System.out.println("Spec1 NumID\t"+specIndexMap1.size());
		System.out.println("Spec2 NumID\t"+specIndexMap2.size());
		System.out.println("Common ID\t"+common.size());
		
		System.out.println("Spec1 only\t"+(specIndexMap1.size()-common.size()));
		for(int specIndex : specIndexMap1.keySet())
		{
			if(common.contains(specIndex))
				continue;
			else
			{
				InsPecTPSM psm = specIndexMap1.get(specIndex);
				System.out.println(specIndex+"\t"+psm.getAnnotation()+"\t"+psm.getCharge()+"\t"+psm.getScore("SpecProb")+"\t"+psm.getScore("FDR"));
			}
		}
		System.out.println("Spec2 only\t"+(specIndexMap2.size()-common.size()));
	}
}
