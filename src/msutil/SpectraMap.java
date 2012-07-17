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
	private ArrayList<Integer> specIndexList = null;
	
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
	
	public synchronized ArrayList<Integer> getSpecIndexList()
	{
		if(specIndexList == null)
		{
			specIndexList = new ArrayList<Integer>(specIndexMap.keySet()); 
			Collections.sort(specIndexList);
		}
		return specIndexList;
	}
}
