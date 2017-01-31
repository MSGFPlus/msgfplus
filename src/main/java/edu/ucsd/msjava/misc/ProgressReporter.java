package edu.ucsd.msjava.misc;

/**
 * @author bryson
 */
public interface ProgressReporter {
        void setProgressData(ProgressData data);
        ProgressData getProgressData();
}
