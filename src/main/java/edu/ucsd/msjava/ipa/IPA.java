package edu.ucsd.msjava.ipa;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Peak;

public class IPA {

    private MS1SpectraMap ms1SpecMap;
    private MSGFPlusResultSet resultSet;
    private Tolerance tol = new Tolerance(5, true);

    public IPA(File deconPeaksFile, File msgfPlusTsvFile) {
        ms1SpecMap = new MS1SpectraMap(deconPeaksFile);
        resultSet = new MSGFPlusResultSet(msgfPlusTsvFile);
    }

    public IPA tolerance(Tolerance tol) {
        this.tol = tol;
        return this;
    }

    public void writeTo(File outputFile) {
        PrintStream out = null;
        try {
            out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

//		class Ion implements Comparable<Ion>
//		{
//			float mz;
//			int charge;
//			float specEValue;
//			
//			@Override
//			public int compareTo(Ion ion) {
//				if(mz > ion.mz)	return 1;
//				else if(mz < ion.mz) return 0;
//				else return charge - ion.charge;
//			}
//			
//			@Override
//			public boolean equals(Object obj)
//			{
//				if(obj instanceof Ion)
//				{
//					Ion ion = (Ion)obj;
//					if(mz == ion.mz && charge == ion.charge)
//						return true;
//				}
//				return false;
//			}
//		}

        HashMap<Integer, HashSet<Float>> precursorMap = new HashMap<Integer, HashSet<Float>>();

        out.println(resultSet.getHeader());
        List<PSM> psmList = resultSet.getPSMList();
        for (PSM psm : psmList) {
            float precursorMz = (float) (psm.getComposition().getMass() / psm.getCharge() + Composition.PROTON);
            int charge = psm.getCharge();

            // 1e-9: Shew, 1e-11: Human

//			if(specEValue < 5e-11 ||  ms1SpecMap.getPrecursorPeaks(psm.getScanNum(), precursorMz, tol) != null 
//			if(eValue < 0.001 || ms1SpecMap.getPrecursorPeaks(psm.getScanNum(), precursorMz, tol) != null 
//					&& ms1SpecMap.getPrecursorPeaks(psm.getScanNum(), secondIsotopeMz, tol) != null
//			if(psm.getScanNum() == 19959)
//				System.out.println("Debug");

//			if(psm.getSpecEValue() < 1e-9 && !ms1SpecMap.checkMS1Peaks(psm.getScanNum(), precursorMz, charge, tol, 3))
//			{
//				System.out.println("\t"+psm.getScanNum()+"\t"+psm.getPeptide() + "\t" + psm.getSpecEValue());
//			}

            if (psm.getSpecEValue() < 1e-10 || ms1SpecMap.checkMS1Peaks(psm.getScanNum(), precursorMz, charge, tol, 3)
                    ) {
                int precursorScan = ms1SpecMap.getPrecursorScan(psm.getScanNum());
                HashSet<Float> identifiedPrecursors = precursorMap.get(precursorScan);
                if (identifiedPrecursors == null) {
                    identifiedPrecursors = new HashSet<Float>();
                    precursorMap.put(precursorScan, identifiedPrecursors);
                }
                if (identifiedPrecursors.contains(precursorMz))
                    continue;
                else
                    identifiedPrecursors.add(precursorMz);

                out.println(psm.getResultString());
            }
        }

        out.close();
    }
}
