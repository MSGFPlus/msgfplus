package edu.ucsd.msjava.ipa;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeMap;

import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Pair;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.BufferedLineReader;

public class MS1SpectraMap {
	private TreeMap<Integer, Spectrum> ms1SpecMap;
	
	public MS1SpectraMap(File peaksFile)
	{
		parsePeaksFile(peaksFile);
	}
	
	public Peak getPrecursorPeaks(int scanNum, float mz, Tolerance tol)
	{
		int precursorScan = ms1SpecMap.lowerKey(scanNum);
		return getMS1Peak(precursorScan, mz, tol);
	}

	public List<Pair<Integer,Float>> getXIC(int scanNum, float mz, Tolerance tol)
	{
		ArrayList<Pair<Integer,Float>> xic = new ArrayList<Pair<Integer,Float>>();
		
		Peak p;
		
		// Move down
		Integer curScanNum = scanNum;
		while(curScanNum != null && scanNum > 0)
		{
			curScanNum = ms1SpecMap.lowerKey(curScanNum);
			if((p = getMS1Peak(curScanNum, mz, tol)) != null)
				xic.add(new Pair<Integer,Float>(curScanNum, p.getIntensity()));
			else
				break;
		}
		
		// Move up
		curScanNum = scanNum;
		while(curScanNum != null && curScanNum < 100000)
		{
			curScanNum = ms1SpecMap.higherKey(curScanNum);
			if((p = getMS1Peak(curScanNum, mz, tol)) != null)
				xic.add(new Pair<Integer,Float>(curScanNum, p.getIntensity()));
			else
				break;
		}
		Collections.sort(xic, new Pair.PairComparator<Integer, Float>());
		return xic;
	}
	
	public Peak getMS1Peak(int ms1ScanNum, float mz, Tolerance tol)
	{
		Spectrum ms1Spec = ms1SpecMap.get(ms1ScanNum);
		if(ms1Spec == null)	return null;
		
		ArrayList<Peak> matchList = ms1Spec.getPeakListByMz(mz, tol);
		if(matchList == null || matchList.size() == 0)
			return null;
		else
		{
			Peak bestPeak = null;
			float distance = Float.MAX_VALUE;
			for(Peak p : matchList)
			{
				float dis = p.getMz() - mz;
				if(dis < 0)
					dis = -dis;
				if(dis < distance)
				{
					distance = dis;
					bestPeak = p;
				}
			}
			return bestPeak;
		}			
	}	
	public boolean checkMS1Peaks(int scanNum, float mz, int charge, Tolerance tol, int windowSize)
	{
		int precursorScanNum = ms1SpecMap.floorKey(scanNum);
		
		
		if(checkMS1Peaks(precursorScanNum, mz, charge, tol))
			return true;
		
		// Move down
		Integer curScanNum = precursorScanNum;
		for(int i=0; i<windowSize; i++)
		{
			curScanNum = ms1SpecMap.lowerKey(curScanNum);
			if(curScanNum == null)
				break;
			else
			{
				if(checkMS1Peaks(curScanNum, mz, charge, tol))
					return true;
			}
		}
		
		// Move up
		curScanNum = precursorScanNum;
		for(int i=0; i<windowSize; i++)
		{
			curScanNum = ms1SpecMap.higherKey(curScanNum);
			if(curScanNum == null)
				break;
			else
			{
				if(checkMS1Peaks(curScanNum, mz, charge, tol))
					return true;
			}
		}
		
		return false;
	}
	
	public boolean checkMS1Peaks(int ms1ScanNum, float mz, int charge, Tolerance tol)
	{
		Spectrum ms1Spec = ms1SpecMap.get(ms1ScanNum);
		if(ms1Spec == null)	return false;
		
		ArrayList<Peak> matchList = ms1Spec.getPeakListByMz(mz, tol);
		if(matchList == null || matchList.size() == 0)
			return false;
		else
		{
			for(Peak p : matchList)
			{
				float secondIsotopeMz = p.getMz() + (float)Composition.ISOTOPE/charge;
				if(ms1Spec.getPeakListByMz(secondIsotopeMz, tol) != null)
					return true;
			}
		}			
		return false;
	}
	
	private void parsePeaksFile(File peaksFile)
	{
		ms1SpecMap = new TreeMap<Integer, Spectrum>();
		
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(peaksFile.getPath());
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		String s;
		String[] token;
		in.readLine();	// header
		int curScanNum = 0;
		Spectrum currentSpec = null;
		while((s=in.readLine()) != null)
		{
			token = s.split("\t");
			int scanNum = Integer.parseInt(token[1]);
			
			if(scanNum > curScanNum)
			{
				if(currentSpec != null)
				{
					ms1SpecMap.put(curScanNum, currentSpec);
				}
				curScanNum = scanNum;
				currentSpec = new Spectrum();
				currentSpec.setMsLevel(1);
			}
			float mz = Float.parseFloat(token[2]);
			float intensity = Float.parseFloat(token[3]);
			currentSpec.add(new Peak(mz, intensity, -1));
		}
	}
}
