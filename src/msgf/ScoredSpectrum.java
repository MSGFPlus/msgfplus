package msgf;

import msutil.Matter;

public interface ScoredSpectrum<T extends Matter> {
	public int getNodeScore(T prm, T srm);
	public int getEdgeScore(T curNode, T prevNode, float edgeMass);
}
