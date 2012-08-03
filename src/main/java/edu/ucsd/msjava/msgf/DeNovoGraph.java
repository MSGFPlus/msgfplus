package edu.ucsd.msjava.msgf;

import java.util.ArrayList;

import edu.ucsd.msjava.msutil.*;

public abstract class DeNovoGraph<T extends Matter> {
	protected T source;
	protected T pmNode;
	protected ArrayList<T> sinkNodes;
	protected ArrayList<T> intermediateNodes;
	
	public T getSource() { return source; }
	public T getPMNode() { return pmNode; }
	public ArrayList<T> getSinkList() { return sinkNodes; }
	public ArrayList<T> getIntermediateNodeList() { return intermediateNodes; }
	
	public abstract boolean isReverse();
	public abstract int getScore(Peptide pep);
	public abstract int getScore(Annotation annotation);
	public abstract int getNodeScore(T node);
//	public abstract int getEdgeScore(T curNode, T prevNode);
	public abstract ArrayList<Edge<T>> getEdges(T curNode);
	public abstract T getComplementNode(T node);
	public abstract AminoAcidSet getAASet();
	
	public static class Edge<T extends Matter> {
		private T prevNode;
		private float probability;
		private int index;	
		private float mass;

		// scores
		private int cleavageScore;
		private int errorScore;
		
		public Edge(T prevNode, float probability, int index, float mass) 
		{
			this.prevNode = prevNode;
			this.probability = probability;
			this.index = index;
			this.mass = mass;
		}
		public T getPrevNode() {
			return prevNode;
		}
		public void setCleavageScore(int cleavageScore)
		{
			this.cleavageScore = cleavageScore;
		}
		public void setErrorScore(int errorScore)
		{
			this.errorScore = errorScore;
		}
		public void setEdgeMass(float mass)
		{
			this.mass = mass;
		}
		public int getEdgeScore() {
			return cleavageScore+errorScore;
		}
		public int getErrorScore() {
			return errorScore;
		}
		public float getEdgeProbability() {
			return probability;
		}
		public int getEdgeIndex() {
			return index;
		}
		public float getEdgeMass() {
			return mass;
		}
	}
}
