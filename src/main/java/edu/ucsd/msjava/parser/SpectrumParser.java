package edu.ucsd.msjava.parser;

import java.util.Map;

import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.SpectrumMetaInfo;

public interface SpectrumParser {
    Spectrum readSpectrum(LineReader lineReader);

    Map<Integer, SpectrumMetaInfo> getSpecMetaInfoMap(BufferedRandomAccessLineReader lineReader);    // specIndex -> filePos
}
