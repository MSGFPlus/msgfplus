package msutil;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import parser.BufferedRandomAccessLineReader;
import parser.SpectrumParser;

public class SpectraMap implements SpectrumAccessorByScanNum {
	private Hashtable<Integer, Long> scanNumMap = null; 	// key: scanNum, value: filePos
	private SpectrumParser parser;
	private BufferedRandomAccessLineReader lineReader;
	private ArrayList<Integer> scanNumList = null;
	
	public SpectraMap(String fileName, SpectrumParser parser)
	{
		lineReader = new BufferedRandomAccessLineReader(fileName);
		
		this.parser = parser;
		// set map
	    scanNumMap = parser.getScanNumMap(lineReader);
	}
	
	public Spectrum getSpectrumByScanNum(int scanNum)
	{
		Long filePos = scanNumMap.get(scanNum);
		if(filePos == null)
			return null;
		else
		{
			lineReader.seek(filePos);
			Spectrum spec = parser.readSpectrum(lineReader);
			spec.setScanNum(scanNum);
			return spec;
		}
	}
	
	public Spectrum getSpectrumByPos(long filePos)
	{
		if(filePos < 0 || filePos > lineReader.size())
			return null;
		else
		{
			lineReader.seek((int)filePos);
			Spectrum spec = parser.readSpectrum(lineReader);
			return spec;
		}
	}
	
	@Override
	public ArrayList<Integer> getScanNumList()
	{
		if(scanNumList == null)
		{
			scanNumList = new ArrayList<Integer>(scanNumMap.keySet()); 
			Collections.sort(scanNumList);
		}
		return scanNumList;
	}
}
