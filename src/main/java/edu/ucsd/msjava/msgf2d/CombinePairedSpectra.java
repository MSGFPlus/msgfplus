package edu.ucsd.msjava.msgf2d;

import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.parser.MgfSpectrumParser;

public class CombinePairedSpectra {
	public static void main(String argv[]) throws Exception
	{
//		addScanNum();
		makePairedSpectra(3);		
		System.out.println("Done");
	}
	
	public static void makePairedSpectra(int charge) throws Exception
	{
		String specFileCID = System.getProperty("user.home")+"/Research/Data/Heck/Tryptic/CID_Scan_03_tryp.mgf";
		String specFileETD = System.getProperty("user.home")+"/Research/Data/Heck/Tryptic/ETD_Scan_03_tryp.mgf";
		
		String specFileOutCID1 = System.getProperty("user.home")+"/Research/MSGF2D/CID_paired_03_tryp_c"+charge+".mgf";
		String specFileOutETD1 = System.getProperty("user.home")+"/Research/MSGF2D/ETD_paired_03_tryp_c"+charge+".mgf";
		
		SpectraMap map1 = new SpectraMap(specFileCID, new MgfSpectrumParser());
		SpectraMap map2 = new SpectraMap(specFileETD, new MgfSpectrumParser());
		
		SpectraContainer container1 = new SpectraContainer();
		SpectraContainer container2 = new SpectraContainer();
		
		for(int scanNum1 : map1.getSpecIndexList())
		{
			int scanNum2 = scanNum1+1;
			Spectrum spec2 = map2.getSpectrumBySpecIndex(scanNum2);
			Spectrum spec1 = map1.getSpectrumBySpecIndex(scanNum1);
			if(spec1 != null && spec2 != null)
			{
				if(spec1.getCharge() == charge && spec2.getCharge() == charge)
				{
					container1.add(spec1);
					container2.add(spec2);
				}
			}
		}
		container1.outputMgfFile(specFileOutCID1);
		container2.outputMgfFile(specFileOutETD1);
	}
	
}
