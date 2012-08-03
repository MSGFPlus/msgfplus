package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Matter;
import edu.ucsd.msjava.msutil.Peak;

public interface ScoredSpectrum<T extends Matter> {
	public int getNodeScore(T prm, T srm);
	public float getNodeScore(T node, boolean isPrefix);
	public int getEdgeScore(T curNode, T prevNode, float edgeMass);
	public boolean getMainIonDirection();	// true: prefix, false: suffix
	public Peak getPrecursorPeak();
	public ActivationMethod[] getActivationMethodArr();
	public int[] getScanNumArr();
}
