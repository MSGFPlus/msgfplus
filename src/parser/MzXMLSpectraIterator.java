package parser;

import java.util.Iterator;

import msutil.Spectrum;


/**
 * A data structure that allows iteration of the mzXML file.
 * Only MS/MS spectra will be returned.
 * @author jung
 *
 */
public class MzXMLSpectraIterator implements Iterator<Spectrum> {
	private MzXMLSpectraMap map;
	private boolean hasNext;
	private int scanNum;
	private Spectrum currentSpectrum;
	
	
	/**
	 * Constructor taking the file name.
	 * @param fileName
	 */
	public MzXMLSpectraIterator(String fileName) {
		this(fileName, 2, Integer.MAX_VALUE);
	}

	/**
	 * Constructor taking the file name and mslevel selectors
	 * @param fileName the path to the mzXML file.
	 * @param minMSLevel spectra with msLevel less than this will be ignored.
	 * @param maxMSLevel spectra with msLevel equal or greater than this will be ignored.
	 */
	public MzXMLSpectraIterator(String fileName, int minMSLevel, int maxMSLevel) {
		map = new MzXMLSpectraMap(fileName).msLevel(minMSLevel, maxMSLevel);
		scanNum = 0;
		currentSpectrum = parseNextSpectrum();
		if(currentSpectrum != null)        hasNext = true;
		else                               hasNext = false;
	}
	
	/**
	 * Get next spectrum.
	 * @return the next spectrum.
	 */
	public Spectrum next() {
		Spectrum curSpecCopy = currentSpectrum;
		currentSpectrum = parseNextSpectrum();
		if(currentSpectrum == null)
			hasNext = false;
		return curSpecCopy;
	}
	
	
	/**
	 * Check whether there is more to parse.
	 * @return true if not done or false 
	 */
	public boolean hasNext() {
		return hasNext;
	}
	
	private Spectrum parseNextSpectrum()
	{
		Spectrum spec = null;
		while(++scanNum <= map.getMaxScanNumber())
		{
			spec = map.getSpectrumBySpecIndex(scanNum);
			if(spec != null)
				break;
		}
		return spec;
	}

	public void remove() {
		assert(false);
	}	
}
