package edu.ucsd.msjava.parser;

import java.util.Collections;
import java.util.Hashtable;
import java.util.Map;

import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.SpectraMap;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.SpectrumMetaInfo;

/**
 * A class that parses MS2 format
 * @author sangtaekim
 *
 */
public class MS2SpectrumParser implements SpectrumParser {

	private Spectrum spec = null;
	private Boolean isSpecSorted = null;
	
	/**
	 * Reads a spectrum from ms2 file and returns it.
	 * @param lineReader A LineReader object points to the start of a spectrum.
	 * @return a spectrum object.
	 */
	public Spectrum readSpectrum(LineReader lineReader) 
	{
		float prevMass = 0;
		String buf;

		do
		{
			buf = lineReader.readLine();
		}
		while(buf != null && buf.startsWith("H"));
		if(buf == null)
			return null;
		
		if(buf.startsWith("S"))
		{
			String[] token = buf.split("\\s+");
			spec = new Spectrum();
			int startScanNum = Integer.parseInt(token[1]);
			int endScanNum = Integer.parseInt(token[2]);
			float precursorMz = Float.parseFloat(token[3]);
			spec = new Spectrum(precursorMz, 0, 0);
			spec.setStartScanNum(startScanNum);
			spec.setEndScanNum(endScanNum);
			isSpecSorted = true;
		}
		else if(spec == null)
			return null;
			
		boolean zParsed = false;
		while((buf = lineReader.readLine()) != null)
		{
			String[] token = buf.split("\\s+");
			if(buf.startsWith("H"))
				continue;
			else if(buf.startsWith("S"))	// start of a next spectrum
			{
				Spectrum specCopy = spec;
				Boolean isSpecSortedCopy = isSpecSorted;
				
				spec = new Spectrum();
				int startScanNum = Integer.parseInt(token[1]);
				int endScanNum = Integer.parseInt(token[2]);
				float precursorMz = Float.parseFloat(token[3]);
				spec = new Spectrum(precursorMz, 0, 0);
				spec.setStartScanNum(startScanNum);
				spec.setEndScanNum(endScanNum);
				isSpecSorted = true;
				
				if(!isSpecSortedCopy)
					Collections.sort(specCopy);
				return specCopy;
            }
			else if(buf.startsWith("Z"))
			{
				if(!zParsed)
				{
					int charge = Integer.parseInt(token[1]);
					float precursorMH = Float.parseFloat(token[2]);
					float precursorMz = ((precursorMH-(float)Composition.PROTON)+charge*(float)Composition.PROTON)/charge;
					spec.setPrecursor(new Peak(precursorMz, 0, charge));
					zParsed = true;
				}
				else
				{
					spec.setPrecursorCharge(0);
				}
			}
			else if(token.length == 2)	// a peak
			{
                assert(spec != null);
                float mass = Float.parseFloat(token[0]);
                if(isSpecSorted && mass < prevMass)
                	isSpecSorted = false;
                else
                	prevMass = mass;
                float intensity = Float.parseFloat(token[1]);
                spec.add(new Peak(mass, intensity, 1));
			}
		}
		
		if(spec != null)
		{
			if(!isSpecSorted)
				Collections.sort(spec);
			Spectrum specCopy = spec;
			spec = null;
			return specCopy;
		}

		return spec;			
	}
	
	/**
	 * Read the entire ms2 file and generates a map from spectrum indexes to file positions of spectra.
	 * @param lineReader A reader points to the start of the spectrum.
	 * @return A Hashtable object maps a spectrum index into a file position.
	 */
	public Map<Integer, SpectrumMetaInfo> getSpecMetaInfoMap(
			BufferedRandomAccessLineReader lineReader) {
		Hashtable<Integer, SpectrumMetaInfo> specIndexMap = new Hashtable<Integer, SpectrumMetaInfo>();
		String buf;
		long offset = 0;
		int specIndex = 0;
		
		SpectrumMetaInfo metaInfo = null;
		while((buf = lineReader.readLine()) != null)
		{
			if(buf.startsWith("S"))	// scan
			{
				specIndex++;
				
				metaInfo = new SpectrumMetaInfo();
				metaInfo.setPosition(offset);
				metaInfo.setID("index=" + (specIndex-1));
				
				String[] token = buf.split("\\s+");
				if(token.length < 4)
				{
					System.err.println("Illegal ms2 file format!");
					System.exit(-1);
				}
				float precursorMz = Float.parseFloat(token[3]);
				metaInfo.setPrecursorMz(precursorMz);
				specIndexMap.put(specIndex, metaInfo);
			}
			
			offset = lineReader.getPosition();
		}
		return specIndexMap;	
	}
	
	public static void test() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/QCShew/QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.ms2";
		
	    java.util.Map<Integer, Float> specIndexPrecursorMzMap = new java.util.HashMap<Integer,Float>();
	    int numSpecs;

	    numSpecs = 0;
	    SpectraMap map = new SpectraMap(fileName, new MS2SpectrumParser());
	    
	    for(int specIndex : map.getSpecIndexList())
	    {
	    	Spectrum spec = map.getSpectrumBySpecIndex(specIndex);
	    	numSpecs++;
	    	specIndexPrecursorMzMap.put(spec.getSpecIndex(), spec.getPrecursorPeak().getMz());
	    }
	    System.out.println("NumSpectra: " + numSpecs);
	    
//	    Spectrum scan87 = map.getSpectrumBySpecIndex(79);
//	    System.out.println("**** " + scan87.getPrecursorPeak().getMz()+" "+scan87.getPrecursorPeak().getCharge());
	    
	    numSpecs = 0;
	    SpectraIterator iterator = new SpectraIterator(fileName, new MS2SpectrumParser());
	    while(iterator.hasNext())
	    {
	    	Spectrum spec = iterator.next();
	    	numSpecs++;
	    	
	    	Float precursorMz = specIndexPrecursorMzMap.get(spec.getSpecIndex());
//	    	System.out.println(spec.getPrecursorPeak().getMz()+" "+spec.getCharge()+" "+spec.getSpecIndex()+" "+spec.getScanNum());
	    	if(precursorMz == null || precursorMz != spec.getPrecursorPeak().getMz())
	    	{
	    		System.out.println(precursorMz + " != " + spec.getPrecursorPeak().getMz());
		    	System.exit(0);
	    	}
	    }

	    System.out.println("NumSpectra: " + numSpecs);
	}
	
	public static void main(String argv[]) throws Exception
	{
		test();
	}
}
