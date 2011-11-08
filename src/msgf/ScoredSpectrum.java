package msgf;

import msutil.ActivationMethod;
import msutil.Matter;
import msutil.Peak;

public interface ScoredSpectrum<T extends Matter> {
	public int getNodeScore(T prm, T srm);
	public float getNodeScore(T node, boolean isPrefix);
	public int getEdgeScore(T curNode, T prevNode, float edgeMass);
	public boolean getMainIonDirection();	// true: prefix, false: suffix
	public Peak getPrecursorPeak();
	public ActivationMethod[] getActivationMethodArr();
	public int[] getScanNumArr();
}
