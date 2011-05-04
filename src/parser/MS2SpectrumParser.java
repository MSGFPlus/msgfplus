package parser;

import java.util.Collections;
import java.util.Hashtable;

import msutil.Composition;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.Spectrum;

/**
 * A class that parses MS2 format
 * @author sangtaekim
 *
 */
public class MS2SpectrumParser implements SpectrumParser {

	/**
	 * Reads a spectrum from ms2 file and returns it.
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
			if(buf.startsWith(":"))	// start of a spectrum
			{
				token = buf.substring(1).split("\\.");
				int startScanNum = Integer.parseInt(token[0]);
				int endScanNum = Integer.parseInt(token[1]);
				
				buf = lineReader.readLine();
				if(buf != null)
				{
					token = buf.split("\\s+");
					float MH = Float.parseFloat(token[0]);
					int charge = Integer.parseInt(token[1]);
					float precursorMz = ((MH-(float)Composition.PROTON)+charge*(float)Composition.PROTON)/charge;
					float precursorIntensity = 0;
	                spec = new Spectrum(precursorMz, charge, precursorIntensity);
	                spec.setStartScanNum(startScanNum);
	                spec.setEndScanNum(endScanNum);
				}
            }
			else if(token.length == 2)	// a peak
			{
                assert(spec != null);
                float mass = Float.parseFloat(token[0]);
                if(sorted && mass < prevMass)
                	sorted = false;
                else
                	prevMass = mass;
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
	 * Read the entire ms2 file and generates a map from scan numbers to file positions of spectra.
	 * @param lineReader A reader points to the start of the spectrum.
	 * @return A Hashtable object maps a scan number into a file position.
	 */
	@Override
	public Hashtable<Integer, Long> getScanNumMap(
			BufferedRandomAccessLineReader lineReader) {
		Hashtable<Integer, Long> scanNumMap = new Hashtable<Integer, Long>();
		String buf;
		long offset = 0;
		while((buf = lineReader.readLine()) != null)
		{
			if(buf.startsWith(":"))
			{
				String[] token = buf.substring(1).split("\\.");
				int startScanNum = Integer.parseInt(token[0]);
                scanNumMap.put(startScanNum, offset);
			}
			
			offset = lineReader.getPosition();
		}
		return scanNumMap;	
	}
	
	public static void test() throws Exception
	{
		// String fileName = System.getProperty("user.home")+"/Research/Data/Shewanella/spectra/SBShew_08m_13Sep04_Andro_0904-1_4-8_dta.ms2";
		String fileName = System.getProperty("user.home")+"/Data/Spectra/SOneSpectra/Yufeng_LTQ_15um_4000psi_50ngShew_dd_4_dta.ms2";
		SpectraIterator iterator = new SpectraIterator(fileName, new MS2SpectrumParser());
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
