package edu.ucsd.msjava.msutil;

import java.util.Collections;


public class TopNFilter implements Reshape {

    // the number of peaks to retain.
    private int topN;


    /**
     * Constructor.
     *
     * @param top the number of peaks to keep.
     */
    public TopNFilter(int topN) {
        this.topN = topN;
    }

    /**
     * Getter method.
     *
     * @return the number of top peaks for this window.
     */
    public int getTop() {
        return topN;
    }


    public Spectrum apply(Spectrum s) {

        Spectrum intSortedSpec = (Spectrum) s.clone();
        Collections.sort(intSortedSpec, Collections.reverseOrder(new Peak.IntensityComparator()));

        Spectrum retSpec = (Spectrum) s.clone();
        retSpec.clear();    // remove all peaks

        for (int peakIndex = 0; peakIndex < topN && peakIndex < intSortedSpec.size(); peakIndex++) {
            retSpec.add(intSortedSpec.get(peakIndex));
        }
        return retSpec;
    }
}
