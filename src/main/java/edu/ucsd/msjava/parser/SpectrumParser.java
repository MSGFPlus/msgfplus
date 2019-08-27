package edu.ucsd.msjava.parser;

import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.SpectrumMetaInfo;

import java.util.Map;

public interface SpectrumParser {

    int MAX_SCAN_MISSING_WARNINGS = 10;

    Spectrum readSpectrum(LineReader lineReader);

    Map<Integer, SpectrumMetaInfo> getSpecMetaInfoMap(BufferedRandomAccessLineReader lineReader);    // specIndex -> filePos

    /**
     * Gets the number of spectra for which the scan number could not be determined
     * @return
     */
    long getScanMissingWarningCount();
}
