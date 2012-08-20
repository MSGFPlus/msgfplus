package edu.ucsd.msjava.msutil;

import java.util.HashMap;
import java.util.Map;

public class SpectrumMetaInfo {
	
	private float precursorMz;
	private String id;
	private long position;	// position in file
	private Map<String,String> additionalMap;
	
	public SpectrumMetaInfo(String id, float precursorMz, long position)
	{
		this.id = id;
		this.precursorMz = precursorMz;
		this.position = position;
	}
	
	public SpectrumMetaInfo() {}
	
	public void setID(String id)	{ this.id = id; }
	public void setPrecursorMz(float precursorMz) { this.precursorMz = precursorMz; }
	public void setPosition(long position) { this.position = position; }
	
	public String getID()
	{
		return id;
	}
	
	public float getPrecursorMz()
	{
		return precursorMz;
	}
	
	public long getPosition()
	{
		return position;
	}
	
	public void setAdditionalInfo(String key, String value)
	{
		if(additionalMap == null)
			additionalMap = new HashMap<String, String>();
		additionalMap.put(key, value);
	}
	
	public String getAdditionalInfo(String key)
	{
		if(additionalMap == null)
			return null;
		else
			return additionalMap.get(key);
	}
}
