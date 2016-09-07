package edu.ucsd.msjava.msscorer;

import edu.ucsd.msjava.msutil.IonType;

public interface NewAdditiveScorer {
    // for scoring nodes
    float getNodeScore(Partition part, IonType ionType, int rank);

    float getMissingIonScore(Partition part, IonType ionType);

    // for scoring edges
    float getErrorScore(Partition part, float error);

    // index => nn:0, ny:1, yn:2, yy:3
    float getIonExistenceScore(Partition part, int index, float probPeak);
}