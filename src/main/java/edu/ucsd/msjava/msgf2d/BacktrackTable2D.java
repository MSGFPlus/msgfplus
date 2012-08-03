package edu.ucsd.msjava.msgf2d;

import java.util.ArrayList;
import java.util.Hashtable;


import edu.ucsd.msjava.msgf.DeNovoGraph;
import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Matter;
import edu.ucsd.msjava.suffixarray.SuffixArray;

public class BacktrackTable2D<T extends Matter> extends Hashtable<T, BacktrackPointer2D> {
	private static final long serialVersionUID = 1L;
	DeNovoGraph<T> graph;
	AminoAcidSet aaSet;

	public BacktrackTable2D(DeNovoGraph<T> graph, AminoAcidSet aaSet)
	{
		this.graph = graph;
		this.aaSet = aaSet;
	}

	public void getReconstructions(T curNode, int score1, int score2, String prefix, ArrayList<String> reconstructions)
	{
		getReconstructions(curNode, score1, score2, prefix, reconstructions, null);
	}

	public void getReconstructions(T curNode, int score1, int score2, String prefix, ArrayList<String> reconstructions, SuffixArray sa)
	{
		//		System.out.println(curNode.getNominalMass()+"\t"+score1+"\t"+score2+"\t"+prefix);
		if(sa != null && sa.search(prefix) < 0)
			return;

		BacktrackPointer2D pointer = this.get(curNode);
		if(pointer == null)
			return;
		if(score1 >= pointer.getMaxScore1() || score2 >= pointer.getMaxScore2())
			return;
		assert(pointer != null);
		if(curNode.equals(graph.getSource()))	// source
		{
			reconstructions.add(prefix);
			return;
		}
		for(DeNovoGraph.Edge<T> edge : graph.getEdges(curNode))
		{
			T prevNode = edge.getPrevNode();
			int edgeIndex = edge.getEdgeIndex();
			if(pointer.isSet(score1, score2, edgeIndex))
				getReconstructions(prevNode, score1-pointer.getCurScore1(), score2-pointer.getCurScore2(), prefix+aaSet.getAminoAcid(edgeIndex).getResidueStr(), reconstructions, sa);
		}
	}

	public String getOneReconstruction(T curNode, int score1, int score2, String prefix)
	{
		BacktrackPointer2D pointer = this.get(curNode);
		if(pointer == null)
			return null;
		if(score1 >= pointer.getMaxScore1() || score2 >= pointer.getMaxScore2())
			return null;
		assert(pointer != null);
		if(curNode.equals(graph.getSource()))	// source
		{
			return prefix;
		}
		for(DeNovoGraph.Edge<T> edge : graph.getEdges(curNode))
		{
			T prevNode = edge.getPrevNode();
			int edgeIndex = edge.getEdgeIndex();
			if(pointer.isSet(score1, score2, edgeIndex))
				getOneReconstruction(prevNode, score1-pointer.getCurScore1(), score2-pointer.getCurScore2(), prefix+aaSet.getAminoAcid(edgeIndex).getResidueStr());
		}
		return null;
	}	
}
