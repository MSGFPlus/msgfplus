package edu.ucsd.msjava.msutil;

public class SpectrumMetaInfo {
	
	private float precursorMz;
	private String id;
	private long position;	// position in file
	
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
}
