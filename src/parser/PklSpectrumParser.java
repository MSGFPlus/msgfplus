package parser;

import java.util.Collections;
import java.util.Hashtable;

import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;

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
		int scanNum = -1;

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
                spec.setScanNum(++scanNum);
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
	 * Read the entire pkl file and generates a map from scan numbers to file positions of spectra.
	 * @param lineReader A reader points to the start of the spectrum.
	 * @return A Hashtable object maps a scan number into a file position.
	 */
	@Override
	public Hashtable<Integer, Long> getScanNumMap(
			BufferedRandomAccessLineReader lineReader) {
		Hashtable<Integer, Long> scanNumMap = new Hashtable<Integer, Long>();
		String buf;
		long offset = 0;
		int sequentialScanNum = -1;
		while((buf = lineReader.readLine()) != null)
		{
			String[] token = buf.split("\\s+");
			if(token.length == 3)	// start of a spectrum
                scanNumMap.put(++sequentialScanNum, offset);
			
			offset = lineReader.getPosition();
		}
		return scanNumMap;	
	}
	
	public static void test() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/NitinSignalPep/HCGS_QToF_Trypsin.pkl";
	    SpectraIterator iterator = new SpectraIterator(fileName, new PklSpectrumParser());
	    int numSpecs = 0;
	    while(iterator.hasNext())
	    {
	    	iterator.next();
	    	numSpecs++;
	    }
	    System.out.println("NumSpectra: " + numSpecs);
	}
	
	public static void main(String argv[]) throws Exception
	{
		test();
	}
}
