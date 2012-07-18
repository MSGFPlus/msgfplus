package jmzparser;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import msutil.SpecFileFormat;

import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.CvParam;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.ParamGroup;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.UserParam;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.ms2_parser.Ms2File;
import uk.ac.ebi.pride.tools.mzdata_parser.MzDataFile;
import uk.ac.ebi.pride.tools.mzml_wrapper.MzMlWrapper;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import uk.ac.ebi.pride.tools.pkl_parser.PklFile;

// Not active because currently jmzreader does not retrive activation method information.
public class SpectraIterator implements Iterator<msutil.Spectrum>, Iterable<msutil.Spectrum> {
	private JMzReader parser;
	private Iterator<uk.ac.ebi.pride.tools.jmzreader.model.Spectrum> itr;
	private int specIndex = 0;
	private SpecFileFormat specFormat;
	
	public SpectraIterator(File specFile, SpecFileFormat specFormat) throws FileNotFoundException
	{
		try {
			if(specFormat == SpecFileFormat.MZML)
			{
				parser = new MzMlWrapper(specFile);
			}
			else if(specFormat == SpecFileFormat.MZXML)
			{
				parser = new MzXMLFile(specFile);
			}
			else if(specFormat == SpecFileFormat.MGF)
			{
				parser = new MgfFile(specFile);
			}
			else if(specFormat == SpecFileFormat.MS2)
			{
				parser = new Ms2File(specFile);
			}
			else if(specFormat == SpecFileFormat.PKL)
			{
				parser = new PklFile(specFile);
			}
			else if(specFormat == SpecFileFormat.MZDATA)
			{
				parser = new MzDataFile(specFile);
			}
			else
			{
				parser = null;
				System.err.println("Unsupported spectrum file format: " + specFile.getPath());
				System.exit(-1);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		itr = parser.getSpectrumIterator();
		this.specFormat = specFormat;
	}
	
	public boolean hasNext() 
	{
		return itr.hasNext();
	}

	public msutil.Spectrum next() 
	{
		uk.ac.ebi.pride.tools.jmzreader.model.Spectrum jmzSpec = itr.next();
		msutil.Spectrum spec = JmzSpectrumConverter.convert(jmzSpec);
		spec.setSpecIndex(++specIndex);
		
		return spec;
	}

	public void remove() 
	{
		throw new UnsupportedOperationException("SpectraIterator.remove() not implemented");
	}
	
	public Iterator<msutil.Spectrum> iterator() {
		return this;
	}
	
	public static void main(String argv[]) throws Exception
	{
		Logger.getRootLogger().setLevel(Level.FATAL);
		test();
	}
	
	public static void test() throws Exception
	{
		File specFile = new File("/Users/kims336/Research/Data/JMzReader/example.mzML");
		specFile = new File("/Users/kims336/Research/Data/JMzReader/small.pwiz.1.1.mzML");

		JMzReader parser = new MzMlWrapper(specFile);
		
//		File specFile = new File("/Users/kims336/Research/Data/JMzReader/example.mzXML");
//		JMzReader parser = new MzXMLFile(specFile);
		
//		JMzReader parser = new MzMlWrapper(specFile);
		
//		File specFile = new File("/Users/kims336/Research/Data/JMzReader/test.mgf");
//		JMzReader parser = new MgfFile(specFile);
		
		Iterator<uk.ac.ebi.pride.tools.jmzreader.model.Spectrum> itr = parser.getSpectrumIterator();
		while(itr.hasNext())
		{
			uk.ac.ebi.pride.tools.jmzreader.model.Spectrum spec = itr.next();
			if(spec.getMsLevel() == 2)
			{
				System.out.println(spec.getId()+" "+spec.getMsLevel()+" "+spec.getPrecursorMZ()+" "+spec.getPrecursorIntensity()+" "+spec.getPrecursorCharge());
				ParamGroup params = spec.getAdditional();
				System.out.println("CVParams");
				for(CvParam param : params.getCvParams())
				{
					System.out.println(param.getName()+"\t|\t"+param.getValue());
				}
				System.out.println("UserParams");
				for(UserParam param : params.getUserParams())
				{
					System.out.println(param.getName()+"\t|\t"+param.getValue());
				}
				System.exit(1);
			}
		}
	}
}