package msgf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;

import suffixarray.SuffixArray;

import msutil.AminoAcidSet;
import msutil.Matter;

public class BacktrackTable<T extends Matter> extends HashMap<T, BacktrackPointer> {
	private static final long serialVersionUID = 1L;
	DeNovoGraph<T> graph;
	
	public BacktrackTable(DeNovoGraph<T> graph)
	{
		this.graph = graph;
	}
	
	public void getReconstructions(T curNode, int score, String prefix, ArrayList<String> reconstructions)
	{
		getReconstructions(curNode, score, prefix, reconstructions, null);
	}
	
	public void getReconstructions(T curNode, int score, String prefix, ArrayList<String> reconstructions, SuffixArray sa)
	{
		if(sa != null && sa.search(prefix) < 0)
			return;
		
		BacktrackPointer pointer = this.get(curNode);
		if(pointer == null)
			return;
		if(score >= pointer.getMaxScore())
			return;
		assert(pointer != null);
		if(curNode.equals(graph.getSource()))	// source
		{
			reconstructions.add(prefix);
			return;
		}
		
		for(DeNovoGraph.Edge<T> edge : graph.getEdges(curNode))
		{
			int edgeIndex = edge.getEdgeIndex();
//			String residue;
//			if(edgeIndex >= 0)
//				residue = String.valueOf(graph.getAASet().getAminoAcid(edgeIndex).getResidue());
//			else
//			{
//				if(edgeIndex == -2)
//					residue="K.";
//				else if(edgeIndex == -3)
//					residue = "G.";
//				else
//					residue = "";
//			}
			if(pointer.isSet(score, edgeIndex))
				getReconstructions(edge.getPrevNode(), score-(edge.getEdgeScore()+pointer.getNodeScore()), prefix+graph.getAASet().getAminoAcid(edgeIndex).getResidueStr(), reconstructions, sa);
		}
	}
	
	public String getOneReconstruction(T curNode, int score, String prefix)
	{
		BacktrackPointer pointer = this.get(curNode);
		if(pointer == null)
			return null;
		if(score >= pointer.getMaxScore())
			return null;
		assert(pointer != null);
		if(curNode.equals(graph.getSource()))	// source
		{
			return prefix;
		}
//		for(T prevNode : graph.getPreviousNodes(curNode))
//		{
//			int edgeIndex = graph.getEdgeIndex(curNode, prevNode);
//			if(pointer.isSet(score, edgeIndex))
//				getOneReconstruction(prevNode, score-pointer.getCurScore(), prefix+aaSet.getAminoAcid(edgeIndex).getResidue());
//		}
		for(DeNovoGraph.Edge<T> edge : graph.getEdges(curNode))
		{
			int edgeIndex = edge.getEdgeIndex();
			if(pointer.isSet(score, edgeIndex))
				getOneReconstruction(edge.getPrevNode(), score-(edge.getEdgeScore()+pointer.getNodeScore()), prefix+graph.getAASet().getAminoAcid(edgeIndex).getResidueStr());
		}		
		return null;
	}	
}
