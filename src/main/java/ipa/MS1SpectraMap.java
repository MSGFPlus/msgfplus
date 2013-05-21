package ipa;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;
import java.util.TreeMap;

import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.BufferedLineReader;

public class MS1SpectraMap {
	TreeMap<Integer, Spectrum> ms1SpecMap;
	
	public MS1SpectraMap(File peaksFile)
	{
		parsePeaksFile(peaksFile);
	}
	
	public List<Peak> getMS1Peaks(int scanNum, float mz, Tolerance tol)
	{
		int precursorScan = ms1SpecMap.floorKey(scanNum);
		
		// TODO: implement it
		return null;
		
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
