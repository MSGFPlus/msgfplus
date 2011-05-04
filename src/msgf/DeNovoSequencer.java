package msgf;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;

import msutil.AminoAcidSet;
import msutil.Matter;
import msutil.Peptide;

public class DeNovoSequencer<T extends Matter> {
	ScoredSpectrum<T> scoredSpec;
	DeNovoGraph<T> graph;
	AminoAcidSet aaSet;
	
	ArrayList<String> deNovoStrings;
	ArrayList<Peptide> deNovoPeptides;
	ArrayList<T> optimalNodes;
	int deNovoScore;
	BacktrackTable<T> backtrackTable;

	public DeNovoSequencer(ScoredSpectrum<T> scoredSpec, DeNovoGraph<T> graph, AminoAcidSet aaSet)
	{
		this.scoredSpec = scoredSpec;
		this.graph = graph;
		this.aaSet = aaSet;
	}

	public DeNovoSequencer(ScoredSpectrum<T> scoredSpec, DeNovoGraph<T> graph)
	{
		this(scoredSpec, graph, AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys());
	}
	
	
	public int getDeNovoScore()	{ return deNovoScore; }
	
	public void deNovoSequencing(boolean backtrack, boolean trypticOnly)
	{
		if(backtrack)
		{
			backtrackTable = new BacktrackTable<T>(graph);
			BacktrackPointer sourcePointer = new BacktrackPointer(0, 1, 0);
			sourcePointer.setBacktrack(0, 0);
			backtrackTable.put(graph.getSource(), sourcePointer);
		}
		
		T source = graph.getSource();
		ArrayList<T> intermediateNodeList = graph.getIntermediateNodeList();
		
		Hashtable<T, Integer> table = new Hashtable<T, Integer>();
		table.put(source, 0);
		
		for(int i=1; i<intermediateNodeList.size(); i++)
		{
			T curNode = intermediateNodeList.get(i);
			T srm = curNode;	// SRM
			T prm = graph.getComplementNode(srm);
			int curScore;
			curScore = scoredSpec.getNodeScore(prm, srm);
			int prevMaxScore = Integer.MIN_VALUE;
			
			BacktrackPointer backPointer = null;
			for(DeNovoGraph.Edge<T> edge : graph.getEdges(curNode))
			{
				T prevNode = edge.getPrevNode();
				Integer prevScore = table.get(prevNode);
				if(prevScore != null)
				{
					if(prevScore > prevMaxScore)
					{
						prevMaxScore = prevScore;
						if(backtrack)
						{
							int curDeNovoScore = prevMaxScore + curScore;
							backPointer = new BacktrackPointer(curDeNovoScore, curDeNovoScore+1, curScore);
							BacktrackPointer prevPointer = backtrackTable.get(prevNode);
							backPointer.addBacktrackPointers(prevPointer, edge.getEdgeIndex(), edge.getEdgeScore());
						}
					}
					else if(prevScore == prevMaxScore)
					{
						if(backtrack)
						{
							assert(backPointer != null);
							BacktrackPointer prevPointer = backtrackTable.get(prevNode);
							backPointer.addBacktrackPointers(prevPointer, edge.getEdgeScore(), edge.getEdgeScore());
						}
					}
				}
			}
			if(prevMaxScore >= Integer.MIN_VALUE)
			{
				if(backtrack)
					backtrackTable.put(curNode, backPointer);
				int curDeNovoScore = prevMaxScore+curScore;
				table.put(curNode, curDeNovoScore);
			}
		}
		
		ArrayList<T> sinkNodeList = graph.getSinkList();
		HashSet<T> optimalNodeSet = null;
		this.deNovoScore = Integer.MIN_VALUE;
		for(T curNode : sinkNodeList)
		{
			BacktrackPointer backPointer = null;
			for(DeNovoGraph.Edge<T> edge : graph.getEdges(curNode))
			{
				T prevNode = edge.getPrevNode();
				Integer prevScore = table.get(prevNode);
				if(prevScore != null)
				{
					if(prevScore > deNovoScore)
					{
						deNovoScore = prevScore;
						optimalNodeSet = new HashSet<T>();
						optimalNodeSet.add(curNode);
						if(backtrack)
						{
							backPointer = new BacktrackPointer(deNovoScore, deNovoScore+1, 0);
							BacktrackPointer prevPointer = backtrackTable.get(prevNode);
							backPointer.addBacktrackPointers(prevPointer, edge.getEdgeIndex(), edge.getEdgeScore());
						}
					}
					else if(prevScore == deNovoScore)
					{
						assert(optimalNodeSet != null);
						optimalNodeSet.add(curNode);
						if(backtrack)
						{
							if(backPointer == null)
								backPointer = new BacktrackPointer(deNovoScore, deNovoScore+1, 0);
							BacktrackPointer prevPointer = backtrackTable.get(prevNode);
							backPointer.addBacktrackPointers(prevPointer, edge.getEdgeIndex(), edge.getEdgeScore());
						}
					}
				}
			}
			if(backtrack && backPointer != null)
				backtrackTable.put(curNode, backPointer);
		}
		if(optimalNodeSet != null)
			this.optimalNodes = new ArrayList<T>(optimalNodeSet);
		table = null;
	}
	
	public ArrayList<String> getDeNovoStrings() 
	{
		if(backtrackTable == null)
			return null;

		if(deNovoStrings == null)
		{
			this.deNovoStrings = new ArrayList<String>();
			if(optimalNodes != null)
			{
				for(T targetNode : optimalNodes)
					backtrackTable.getReconstructions(targetNode, deNovoScore, "", deNovoStrings);
			}
		}
		
		return deNovoStrings; 
	}
	
	public ArrayList<Peptide> getDeNovoPeptides()
	{
		if(deNovoPeptides != null)
			return deNovoPeptides;
		else 
		{
			if(deNovoStrings == null)
				getDeNovoStrings();
			if(deNovoStrings != null)
			{
				deNovoPeptides = new ArrayList<Peptide>();
				for(String str : deNovoStrings)
					deNovoPeptides.add(new Peptide(str, aaSet));
				return deNovoPeptides;
			}
		}
		return null;
	}
	
	public ArrayList<T> getOptimalNodes()	{ return optimalNodes; }
}
