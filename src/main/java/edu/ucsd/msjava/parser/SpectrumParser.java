package edu.ucsd.msjava.parser;

import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.SpectrumMetaInfo;

import java.util.Map;

public interface SpectrumParser {
    Spectrum readSpectrum(LineReader lineReader);

    Map<Integer, SpectrumMetaInfo> getSpecMetaInfoMap(BufferedRandomAccessLineReader lineReader);    // specIndex -> filePos
}
