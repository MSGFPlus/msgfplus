package jmzparser;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;

import msutil.SpecFileFormat;

import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.ms2_parser.Ms2File;
import uk.ac.ebi.pride.tools.mzdata_parser.MzDataFile;
import uk.ac.ebi.pride.tools.mzml_wrapper.MzMlWrapper;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import uk.ac.ebi.pride.tools.pkl_parser.PklFile;

public class MzMLSpectraIterator implements Iterator<msutil.Spectrum>, Iterable<msutil.Spectrum> {
	private MzMLUnmarshaller unmarshaller;
	private MzMLObjectIterator<uk.ac.ebi.jmzml.model.mzml.Spectrum> itr;
	private int specIndex = 0;
	
	public MzMLSpectraIterator(File specFile) throws FileNotFoundException
	{
		unmarshaller = new MzMLUnmarshaller(specFile);
		itr = unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", uk.ac.ebi.jmzml.model.mzml.Spectrum.class);
	}
	
	@Override
	public boolean hasNext() 
	{
		return itr.hasNext();
	}

	@Override
	public msutil.Spectrum next() 
	{
		uk.ac.ebi.jmzml.model.mzml.Spectrum jmzSpec = itr.next();
		return null;
//		msutil.Spectrum spec = JmzSpectrumConverter.convert(jmzSpec);
//		spec.setSpecIndex(++specIndex);
//		
//		return spec;
	}

	@Override
	public void remove() 
	{
		throw new UnsupportedOperationException("SpectraIterator.remove() not implemented");
	}
	
	@Override
	public Iterator<msutil.Spectrum> iterator() {
		return this;
	}
}