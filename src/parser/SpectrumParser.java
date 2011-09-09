package parser;

import java.util.Hashtable;

import msutil.Spectrum;

public interface SpectrumParser {
	public Spectrum readSpectrum(LineReader lineReader);
	Hashtable<Integer, Long> getSpecIndexMap(BufferedRandomAccessLineReader lineReader);	// specIndex -> filePos
}
