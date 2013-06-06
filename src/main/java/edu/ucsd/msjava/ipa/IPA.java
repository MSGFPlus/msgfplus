package edu.ucsd.msjava.ipa;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;

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

		class Ion implements Comparable<Ion>
		{
			float mz;
			int charge;
			float specEValue;
			
			@Override
			public int compareTo(Ion ion) {
				if(mz > ion.mz)	return 1;
				else if(mz < ion.mz) return 0;
				else return charge - ion.charge;
			}
			
			@Override
			public boolean equals(Object obj)
			{
				if(obj instanceof Ion)
				{
					Ion ion = (Ion)obj;
					if(mz == ion.mz && charge == ion.charge)
						return true;
				}
				return false;
			}
		}
		
		TreeMap<Float,Ion> precursorMap = new TreeMap<Float,Ion>();
		
		out.println(resultSet.getHeader());
		List<PSM> psmList = resultSet.getPSMList();
		for(PSM psm : psmList)
		{
			float precursorMz = (float)(psm.getComposition().getMass()/psm.getCharge() + Composition.H);
			float secondIsotopeMz = (float)((psm.getComposition().getMass()+Composition.ISOTOPE)/psm.getCharge() + Composition.H);
			float specEValue = psm.getSpecEValue();
			float eValue = psm.getEValue();
			// 1e-9: Shew, 1e-11: Human
//			if(specEValue < 1e-11 || ms1SpecMap.getPrecursorPeaks(psm.getScanNum(), precursorMz, tol) != null 
//			if(eValue < 0.001 || ms1SpecMap.getPrecursorPeaks(psm.getScanNum(), precursorMz, tol) != null 
//					&& ms1SpecMap.getPrecursorPeaks(psm.getScanNum(), secondIsotopeMz, tol) != null
//				)
			if(ms1SpecMap.checkMS1Peaks(psm.getScanNum(), precursorMz, psm.getCharge(), tol, 5))
			{
				out.println(psm.getResultString());
			}
		}
		
		out.close();
	}
}
