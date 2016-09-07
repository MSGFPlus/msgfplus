package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.Matter;

import java.util.List;

public class ScoredSpectrumSumHP<T extends Matter> extends ScoredSpectrumSum<T> {
    public ScoredSpectrumSumHP(List<ScoredSpectrum<T>> scoredSpecList) {
        super(scoredSpecList);
    }


}
