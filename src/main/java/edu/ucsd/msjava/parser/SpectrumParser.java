package edu.ucsd.msjava.parser;

import java.util.Map;

import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.SpectrumMetaInfo;

public interface SpectrumParser {
	public Spectrum readSpectrum(LineReader lineReader);
	Map<Integer, SpectrumMetaInfo> getSpecIndexMap(BufferedRandomAccessLineReader lineReader);	// specIndex -> filePos
}
