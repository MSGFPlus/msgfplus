package edu.ucsd.msjava.ipa;

import java.io.File;

public class IPA {
	
	private MS1SpectraMap ms1SpecMap;
	private MSGFPlusResultMap resultMap;
	
	public IPA(File deconPeaksFile, File msgfPlusTsvFile)
	{
		ms1SpecMap = new MS1SpectraMap(deconPeaksFile);
		resultMap = new MSGFPlusResultMap(msgfPlusTsvFile);
	}
	
	public void filterByMS1Mass()
	{
		
	}
}
