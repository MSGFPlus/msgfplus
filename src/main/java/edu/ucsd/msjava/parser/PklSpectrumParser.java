package edu.ucsd.msjava.parser;

import java.util.Collections;
import java.util.Hashtable;
import java.util.Map;

import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.SpectraMap;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.SpectrumMetaInfo;

/**
 * A class that parses Pkl format
 * @author sangtaekim
 *
 */
public class PklSpectrumParser implements SpectrumParser {

	/**
	 * Reads a spectrum from pkl file and returns it.
	 * @param lineReader A LineReader object points to the start of a spectrum.
	 * @return a spectrum object.
	 */
	public Spectrum readSpectrum(LineReader lineReader) 
	{
		Spectrum spec = null;

		boolean sorted = true;
		float prevMass = 0;
		
		String buf;
		while((buf = lineReader.readLine()) != null)
		{
			String[] token = buf.split("\\s+");
			if(token.length == 3)	// start of a spectrum
			{
				float precursorMz = Float.parseFloat(token[0]);
				float precursorIntensity = Float.parseFloat(token[1]);
				int charge = Integer.parseInt(token[2]);
                spec = new Spectrum(precursorMz, charge, precursorIntensity);
            }
			else if(token.length == 2)	// a peak
			{
                assert(spec != null);
                float mass = Float.parseFloat(token[0]);
                if(sorted && mass < prevMass)
                	sorted = false;
                else
                	prevMass = mass;
//                if(token[1].endsWith("null"))
//                	token[1] = token[1].substring(0, token[1].lastIndexOf("null"));
                float intensity = Float.parseFloat(token[1]);
                spec.add(new Peak(mass, intensity, 1));
			}
  			else 	// end of a spectrum
  			{
  				if(spec != null)
  				{
  	  				if(!sorted)
  	  					Collections.sort(spec);
  					return spec;
  				}
  			}
		}
		return spec;			
	}
	
	/**
	 * Read the entire pkl file and generates a map from spectrum indexes to file positions of spectra.
	 * @param lineReader A reader points to the start of the spectrum.
	 * @return A Hashtable object maps a spectrum index into a file position.
	 */
	public Map<Integer, SpectrumMetaInfo> getSpecIndexMap(
			BufferedRandomAccessLineReader lineReader) {
		Hashtable<Integer, SpectrumMetaInfo> specIndexMap = new Hashtable<Integer, SpectrumMetaInfo>();
		String buf;
		long offset = 0;
		int specIndex = 0;
		while((buf = lineReader.readLine()) != null)
		{
			String[] token = buf.split("\\s+");
			if(token.length == 3)	// start of a spectrum
			{
//                specIndexMap.put(++specIndex, offset);
				++specIndex;
				float precursorMz = Float.parseFloat(token[0]);
				SpectrumMetaInfo metaInfo = new SpectrumMetaInfo();
				metaInfo.setID("index="+(specIndex-1));
				metaInfo.setPrecursorMz(precursorMz);
				metaInfo.setPosition(offset);
				specIndexMap.put(specIndex, metaInfo);
			}
			
			offset = lineReader.getPosition();
		}
		return specIndexMap;	
	}
	
	public static void test() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/ToolDistribution/RefTest/SpecFormatTest/TestSpectra.pkl";
	    SpectraIterator iterator = new SpectraIterator(fileName, new PklSpectrumParser());
	    int numSpecs = 0;
	    while(iterator.hasNext())
	    {
	    	Spectrum spec = iterator.next();
	    	numSpecs++;
	    	System.out.println(spec.getPrecursorPeak().getMz()+" "+spec.getCharge()+" "+spec.getSpecIndex()+" "+spec.getScanNum());
	    }
	    System.out.println("NumSpectra: " + numSpecs);
	    
	    numSpecs = 0;
	    SpectraMap map = new SpectraMap(fileName, new PklSpectrumParser());
	    for(int specIndex : map.getSpecIndexList())
	    {
	    	Spectrum spec = map.getSpectrumBySpecIndex(specIndex);
	    	numSpecs++;
	    	System.out.println(spec.getPrecursorPeak().getMz()+" "+spec.getCharge()+" "+spec.getSpecIndex()+" "+spec.getScanNum());
	    }
	    System.out.println("NumSpectra: " + numSpecs);
	}
	
	public static void main(String argv[]) throws Exception
	{
		test();
	}
}
