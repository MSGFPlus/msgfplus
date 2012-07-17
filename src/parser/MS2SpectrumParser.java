package parser;

import java.util.Collections;
import java.util.Hashtable;

import msutil.Composition;
import msutil.Peak;
import msutil.SpectraIterator;
import msutil.SpectraMap;
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
	 * Read the entire ms2 file and generates a map from spectrum indexes to file positions of spectra.
	 * @param lineReader A reader points to the start of the spectrum.
	 * @return A Hashtable object maps a spectrum index into a file position.
	 */
	public Hashtable<Integer, Long> getSpecIndexMap(
			BufferedRandomAccessLineReader lineReader) {
		Hashtable<Integer, Long> specIndexMap = new Hashtable<Integer, Long>();
		String buf;
		long offset = 0;
		int specIndex = 0;
		while((buf = lineReader.readLine()) != null)
		{
			if(buf.startsWith(":"))
			{
				specIndex++;
				specIndexMap.put(specIndex, offset);
//				String[] token = buf.substring(1).split("\\.");
//				int startScanNum = Integer.parseInt(token[0]);
//                scanNumMap.put(startScanNum, offset);
			}
			
			offset = lineReader.getPosition();
		}
		return specIndexMap;	
	}
	
	public static void test() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/ToolDistribution/RefTest/SpecFormatTest/Yufeng_8LTQ005_010_F3a_dta.ms2";
	    SpectraIterator iterator = new SpectraIterator(fileName, new MS2SpectrumParser());
	    int numSpecs = 0;
	    while(iterator.hasNext())
	    {
	    	Spectrum spec = iterator.next();
	    	numSpecs++;
	    	System.out.println(spec.getPrecursorPeak().getMz()+" "+spec.getCharge()+" "+spec.getSpecIndex()+" "+spec.getScanNum());
	    	if(numSpecs > 5)
	    		break;
	    }
	    System.out.println("NumSpectra: " + numSpecs);
	    
	    numSpecs = 0;
	    SpectraMap map = new SpectraMap(fileName, new MS2SpectrumParser());
	    for(int specIndex : map.getSpecIndexList())
	    {
	    	Spectrum spec = map.getSpectrumBySpecIndex(specIndex);
	    	numSpecs++;
	    	System.out.println(spec.getPrecursorPeak().getMz()+" "+spec.getCharge()+" "+spec.getSpecIndex()+" "+spec.getScanNum());
	    	if(numSpecs > 5)
	    		break;
	    }
	    System.out.println("NumSpectra: " + numSpecs);
	}
	
	public static void main(String argv[]) throws Exception
	{
		test();
	}
}
