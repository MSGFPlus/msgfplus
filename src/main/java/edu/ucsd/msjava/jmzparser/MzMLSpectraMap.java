package edu.ucsd.msjava.jmzparser;

import java.util.ArrayList;

import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshallerException;

import edu.ucsd.msjava.msutil.SpectrumAccessorBySpecIndex;

public class MzMLSpectraMap implements SpectrumAccessorBySpecIndex 
{
	private final MzMLUnmarshaller unmarshaller;
	private final int minMSLevel;		// inclusive
	private final int maxMSLevel;		// exclusive
	
	public MzMLSpectraMap(MzMLAdapter mzmlAdapter)
	{
		unmarshaller = mzmlAdapter.getUnmarshaller();
		minMSLevel = mzmlAdapter.getMinMSLevel();
		maxMSLevel = mzmlAdapter.getMaxMSLevel();
	}

	public edu.ucsd.msjava.msutil.Spectrum getSpectrumBySpecIndex(int specIndex) 
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
				edu.ucsd.msjava.msutil.Spectrum spec = SpectrumConverter.getSpectrumFromJMzMLSpec(jmzSpec);
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

	@Override
	public String getID(int specIndex) {
		return unmarshaller.getSpectrumIDFromSpectrumIndex(specIndex-1);	
	}

	@Override
	public Float getPrecursorMz(int specIndex) {
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
				float precursorMz = SpectrumConverter.getPrecursorMzFromJMzMLSpec(jmzSpec);
				return precursorMz;
			}
		}
		return null;
	}

	@Override
	public String getTitle(int specIndex) {
		return null;
	}
}