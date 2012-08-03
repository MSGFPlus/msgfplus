package edu.ucsd.msjava.parser;

import java.util.Hashtable;

import edu.ucsd.msjava.msutil.Spectrum;

public interface SpectrumParser {
	public Spectrum readSpectrum(LineReader lineReader);
	Hashtable<Integer, Long> getSpecIndexMap(BufferedRandomAccessLineReader lineReader);	// specIndex -> filePos
}
