package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.Matter;

public class ProfilePeak<T extends Matter> implements Comparable<ProfilePeak<T>> {
	T node;
	float probability;
	
	public ProfilePeak(T node, float probability) {
		this.node = node;
		this.probability = probability;
	}
	public T getNode() {
		return node;
	}
	public void setNode(T node) {
		this.node = node;
	}
	public float getProbability() {
		return probability;
	}
	public void setProbability(float probability) {
		this.probability = probability;
	}
	public int compareTo(ProfilePeak<T> p) {
		return node.compareTo(p.node);
	}
}