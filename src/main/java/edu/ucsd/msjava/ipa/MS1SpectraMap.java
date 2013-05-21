package edu.ucsd.msjava.ipa;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Map.Entry;
import java.util.TreeMap;

import edu.ucsd.msjava.msgf.Tolerance;
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
		Entry<Integer, Spectrum> precursorEntry = ms1SpecMap.floorEntry(scanNum);
		
//		int precursorScan = precursorEntry.getKey();
		Spectrum ms1Spec = precursorEntry.getValue();
		if(ms1Spec == null)
		{
			System.out.println("MS1 spec null: " + scanNum);
			return null;
		}
		
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
