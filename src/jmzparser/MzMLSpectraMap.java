package jmzparser;

import java.io.File;
import java.util.ArrayList;

import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshallerException;

import msutil.SpectrumAccessorBySpecIndex;

public class MzMLSpectraMap implements SpectrumAccessorBySpecIndex 
{
	private MzMLUnmarshaller unmarshaller;
	private int minMSLevel = 2;		// inclusive
	private int maxMSLevel = Integer.MAX_VALUE;		// exclusive
	
	public MzMLSpectraMap(File specFile)
	{
		MzMLSpectraIterator.turnOffLogs();
		unmarshaller = new MzMLUnmarshaller(specFile);
	}

	public MzMLSpectraMap(String specFileName)
	{
		unmarshaller = new MzMLUnmarshaller(new File(specFileName));
	}
	
	/**
	 * Setter to set msLevel.
	 * @param minMSLevel minimum msLevel to be considered (inclusive).
	 * @param maxMSLevel maximum msLevel to be considered (inclusive).
	 * @return this object.
	 */
	public MzMLSpectraMap msLevel(int minMSLevel, int maxMSLevel) 
	{ 
		this.minMSLevel = minMSLevel; 
		this.maxMSLevel = maxMSLevel; 
		return this; 
	}
	
	public msutil.Spectrum getSpectrumBySpecIndex(int specIndex) 
	{
		String specID = unmarshaller.getSpectrumIDFromSpectrumIndex(specIndex-1);
		if(specID != null)
		{
			uk.ac.ebi.jmzml.model.mzml.Spectrum jmzSpec = null;
			try {
				jmzSpec = unmarshaller.getSpectrumById(specID);
			} catch (MzMLUnmarshallerException e) {
				e.printStackTrace();
			}
			if(jmzSpec != null)
			{
				msutil.Spectrum spec = MzMLSpectraIterator.getSpectrumFromJMzMLSpec(jmzSpec);
				if(spec.getMSLevel() < minMSLevel || spec.getMSLevel() > maxMSLevel)
					return null;
				else
					return spec; 
			}
		}
		return null;
	}

	public ArrayList<Integer> getSpecIndexList() 
	{
		return new ArrayList<Integer>(unmarshaller.getSpectrumIndexes());
	}
}