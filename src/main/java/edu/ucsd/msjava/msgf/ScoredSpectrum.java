package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Matter;
import edu.ucsd.msjava.msutil.Peak;

public interface ScoredSpectrum<T extends Matter> {
    int getNodeScore(T prm, T srm);

    float getNodeScore(T node, boolean isPrefix);

    int getEdgeScore(T curNode, T prevNode, float edgeMass);

    boolean getMainIonDirection();    // true: prefix, false: suffix

    Peak getPrecursorPeak();

    ActivationMethod[] getActivationMethodArr();

    int[] getScanNumArr();
}
