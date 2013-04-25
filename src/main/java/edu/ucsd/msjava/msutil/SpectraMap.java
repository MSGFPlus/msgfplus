package edu.ucsd.msjava.msutil;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;

import edu.ucsd.msjava.parser.BufferedRandomAccessLineReader;
import edu.ucsd.msjava.parser.SpectrumParser;

public class SpectraMap implements SpectrumAccessorBySpecIndex {
	private Map<Integer, SpectrumMetaInfo> specIndexMap = null; 	// key: specIndex, value: metaInfo
	private SpectrumParser parser;
	protected BufferedRandomAccessLineReader lineReader;
	private ArrayList<Integer> specIndexList = null;
	
	private Map<String, Integer> idToIndex = null;
	
	public SpectraMap(String fileName, SpectrumParser parser)
	{
		lineReader = new BufferedRandomAccessLineReader(fileName);
		
		this.parser = parser;
		// set map
	    specIndexMap = parser.getSpecMetaInfoMap(lineReader);
	}
	
	@Override
	public Spectrum getSpectrumById(String specId)
	{
		if(idToIndex == null)
			makeIdToIndexMap();
		Integer specIndex = idToIndex.get(specId);
		if(specIndex == null)
			return null;
		else
			return getSpectrumBySpecIndex(specIndex);
	}
	
	@Override
	public synchronized Spectrum getSpectrumBySpecIndex(int specIndex)
	{
		Long filePos = getFileOffset(specIndex);
		if(filePos == null)
			return null;
		else
		{
			lineReader.seek(filePos);
			Spectrum spec = parser.readSpectrum(lineReader);
			spec.setSpecIndex(specIndex);
			spec.determineIsCentroided();
			spec.setID("index="+String.valueOf(specIndex-1));
			return spec;
		}
	}

	@Override
	public Float getPrecursorMz(int specIndex)
	{
		SpectrumMetaInfo metaInfo = specIndexMap.get(specIndex);
		if(metaInfo == null)
			return null;
		else
			return metaInfo.getPrecursorMz();
	}
	
	@Override
	public String getID(int specIndex)
	{
		SpectrumMetaInfo metaInfo = specIndexMap.get(specIndex);
		if(metaInfo == null)
			return null;
		else
			return metaInfo.getID();
	}

	@Override
	public String getTitle(int specIndex) {
		SpectrumMetaInfo metaInfo = specIndexMap.get(specIndex);
		if(metaInfo == null)
			return null;
		else
			return metaInfo.getAdditionalInfo("title");
	}
	
	public Long getFileOffset(int specIndex)
	{
		SpectrumMetaInfo metaInfo = specIndexMap.get(specIndex);
		if(metaInfo == null)
			return null;
		else
			return metaInfo.getPosition();
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
	
	private void makeIdToIndexMap()
	{
		idToIndex = new HashMap<String, Integer>();
		Iterator<Entry<Integer, SpectrumMetaInfo>> itr = specIndexMap.entrySet().iterator();
		while(itr.hasNext())
		{
			Entry<Integer, SpectrumMetaInfo> entry = itr.next();
			idToIndex.put(entry.getValue().getID(), entry.getKey());
		}
	}
	
}
