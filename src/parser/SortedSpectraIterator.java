package parser;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import msutil.Pair;
import msutil.Spectrum;
import msutil.SpectrumAccessorBySpecIndex;


/**
 * A data structure that allows iteration of the spectrum file by increasing order of parent mass.
 * @author sangtae
 *
 */
public class SortedSpectraIterator implements Iterator<Spectrum> {
	private SpectrumAccessorBySpecIndex map;
	private boolean hasNext;
	private Spectrum currentSpectrum;
	private ArrayList<Integer> specIndexList;
	private int index;
	private final int numSpecs;
	
	/**
	 * Constructor taking the file name.
	 * @param fileName
	 */
	public SortedSpectraIterator(Iterator<Spectrum> itr, SpectrumAccessorBySpecIndex map) {
		int numSpecs = 0;
		ArrayList<Pair<Integer,Float>> scanNumPMPairList = new ArrayList<Pair<Integer,Float>>();
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			scanNumPMPairList.add(new Pair<Integer,Float>(spec.getSpecIndex(), spec.getParentMass()));
			numSpecs++;
		}
		this.numSpecs = numSpecs;
		Collections.sort(scanNumPMPairList, new Pair.PairComparator<Integer,Float>(true));
		specIndexList = new ArrayList<Integer>();
		for(Pair<Integer,Float> p : scanNumPMPairList)
			specIndexList.add(p.getFirst());
		
		this.map = map;
		index = -1;
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
	
	/**
	 * Returns the number of spectra.
	 * @return the number of spectra.
	 */
	public int size() {
		return numSpecs;
	}
	
	private Spectrum parseNextSpectrum()
	{
		++index;
		if(index >= specIndexList.size())
			return null;
		int specIndex = specIndexList.get(index);
		return map.getSpectrumBySpecIndex(specIndex);
	}

	public void remove() {
		assert(false);
	}	
	
	public static void main(String argv[])
	{
		long time = System.currentTimeMillis();
		String fileName = "/home/sangtaekim/Research/Data/HeckWhole/Spectra/090121_NM_Trypsin_20.mzXML";
		MzXMLSpectraIterator itr = new MzXMLSpectraIterator(fileName);
		MzXMLSpectraMap map = new MzXMLSpectraMap(fileName);
		SortedSpectraIterator sortedItr = new SortedSpectraIterator(itr, map);
		int numSpecs = 0;
		while(sortedItr.hasNext())
		{
			Spectrum spec = sortedItr.next();
			System.out.println(spec.getParentMass());
			numSpecs++;
		}
		System.out.println("NumSpecs: " + numSpecs);
		System.out.println("Time: " + (System.currentTimeMillis()-time));
	}
}
