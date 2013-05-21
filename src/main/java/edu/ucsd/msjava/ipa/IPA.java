package edu.ucsd.msjava.ipa;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.List;

import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.Composition;

public class IPA {
	
	private MS1SpectraMap ms1SpecMap;
	private MSGFPlusResultSet resultSet;
	private Tolerance tol = new Tolerance(5, true);
	
	public IPA(File deconPeaksFile, File msgfPlusTsvFile)
	{
		ms1SpecMap = new MS1SpectraMap(deconPeaksFile);
		resultSet = new MSGFPlusResultSet(msgfPlusTsvFile);
	}
	
	public IPA tolerance(Tolerance tol)
	{
		this.tol = tol;
		return this;
	}
	
	public void writeTo(File outputFile)
	{
		PrintStream out = null;
		try {
			out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		out.println(resultSet.getHeader());
		List<PSM> psmList = resultSet.getPSMList();
		for(PSM psm : psmList)
		{
			float precursorMz = (float)(psm.getComposition().getMass()/psm.getCharge() + Composition.H);
			float secondIsotopeMz = (float)((psm.getComposition().getMass()+Composition.ISOTOPE)/psm.getCharge() + Composition.H);
			float specEValue = psm.getSpecEValue();
			if(specEValue < 1e-10 || ms1SpecMap.getPrecursorPeaks(psm.getScanNum(), precursorMz, tol) != null 
//					&& ms1SpecMap.getPrecursorPeaks(psm.getScanNum(), secondIsotopeMz, tol) != null
				)
			{
				out.println(psm.getResultString());
			}
		}
		
		out.close();
	}
}
