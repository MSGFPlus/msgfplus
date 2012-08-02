package jmzparser;

import java.io.File;
import java.util.ArrayList;

import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshallerException;

import msutil.SpectrumAccessorBySpecIndex;

public class MzMLSpectraMap implements SpectrumAccessorBySpecIndex 
{
	private MzMLUnmarshaller unmarshaller;
	
	public MzMLSpectraMap(File specFile)
	{
		unmarshaller = new MzMLUnmarshaller(specFile);
	}

	public MzMLSpectraMap(String specFileName)
	{
		unmarshaller = new MzMLUnmarshaller(new File(specFileName));
	}
	
	public msutil.Spectrum getSpectrumBySpecIndex(int specIndex) 
	{
		
		String specID = unmarshaller.getSpectrumIDFromSpectrumIndex(specIndex);
		if(specID != null)
		{
			uk.ac.ebi.jmzml.model.mzml.Spectrum jmzSpec = null;
			try {
				jmzSpec = unmarshaller.getSpectrumById(specID);
			} catch (MzMLUnmarshallerException e) {
				e.printStackTrace();
			}
			if(jmzSpec != null)
				return MzMLSpectraIterator.getSpectrumFromJMzMLSpec(jmzSpec);
		}
		return null;
	}

	public ArrayList<Integer> getSpecIndexList() 
	{
		return new ArrayList<Integer>(unmarshaller.getSpectrumIndexes());
	}
}