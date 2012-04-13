package msutil;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import parser.BufferedRandomAccessLineReader;
import parser.SpectrumParser;

public class SpectraMap implements SpectrumAccessorBySpecIndex {
	private Hashtable<Integer, Long> specIndexMap = null; 	// key: specIndex, value: filePos
	private SpectrumParser parser;
	protected BufferedRandomAccessLineReader lineReader;
	private ArrayList<Integer> scanNumList = null;
	
	public SpectraMap(String fileName, SpectrumParser parser)
	{
		lineReader = new BufferedRandomAccessLineReader(fileName);
		
		this.parser = parser;
		// set map
	    specIndexMap = parser.getSpecIndexMap(lineReader);
	}
	
	public synchronized Spectrum getSpectrumBySpecIndex(int specIndex)
	{
		Long filePos = specIndexMap.get(specIndex);
		if(filePos == null)
			return null;
		else
		{
			lineReader.seek(filePos);
			Spectrum spec = parser.readSpectrum(lineReader);
			spec.setSpecIndex(specIndex);
			return spec;
		}
	}
	
	@Override
	public synchronized ArrayList<Integer> getSpecIndexList()
	{
		if(scanNumList == null)
		{
			scanNumList = new ArrayList<Integer>(specIndexMap.keySet()); 
			Collections.sort(scanNumList);
		}
		return scanNumList;
	}
}
