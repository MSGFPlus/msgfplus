package msgf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Enzyme;
import msutil.Matter;
import msutil.Peptide;
import msutil.Sequence;

public abstract class MassFactory<T extends Matter> implements DeNovoNodeFactory<T> {

	protected AminoAcidSet aaSet;
	protected ArrayList<T> allNodes;
	protected HashMap<T, ArrayList<DeNovoGraph.Edge<T>>> edgeMap;
	protected Enzyme enzyme;
	protected int maxLength;
	
	public MassFactory(AminoAcidSet aaSet, Enzyme enzyme, int maxLength)
	{
		this.aaSet = aaSet;
		this.enzyme = enzyme;
		this.maxLength = maxLength;
	}

	// true if this graph represents reverse peptides
	@Override
	public boolean isReverse()
	{
		return enzyme == null || enzyme.isCTerm();
	}
	
	public int getMaxLength()	{ return maxLength; }
	
	public ArrayList<T> getAllNodes()	{ return allNodes; }
	
	@Override
	public int size() {
		return allNodes.size(); 
	}
	
	@Override
	public AminoAcidSet getAASet() {
		return aaSet;
	}

//	@Override
//	public ArrayList<DeNovoGraph.Edge<T>> getEdges(T curNode)
//	{
//		return edgeMap.get(curNode);
//	}
	
	@Override
	public DeNovoGraph.Edge<T> getEdge(T curNode, T prevNode)
	{
		for(DeNovoGraph.Edge<T> edge : getEdges(curNode))
		{
			if(edge.getPrevNode().equals(prevNode))
				return edge;
		}
		return null;
	}
	
	@Override
	public T getPreviousNode(T curNode, AminoAcid aa) 
	{
		int aaIndex = aaSet.getIndex(aa);
		for(DeNovoGraph.Edge<T> edge : getEdges(curNode))
		{
			if(edge.getEdgeIndex() == aaIndex)
				return edge.getPrevNode();
		}
		return null;
	}
	
	@Override
	public abstract T getZero();
	
	@Override
	public Enzyme getEnzyme()
	{
		return enzyme;
	}
	
	@Override
	public ArrayList<T> getIntermediateNodes(ArrayList<T> destNodes) 
	{
		HashSet<T> effectiveNodeSet = new HashSet<T>();
		ArrayList<T> curFreshNodes = new ArrayList<T>(destNodes);
		
		while(!curFreshNodes.isEmpty())
		{
			ArrayList<T> newFreshNodes = new ArrayList<T>();
			for(T node : curFreshNodes)
			{
				ArrayList<DeNovoGraph.Edge<T>> edges = getEdges(node);
				if(edges != null)
				{
					for(DeNovoGraph.Edge<T> edge : edges)
					{
						T prevNode = edge.getPrevNode();
						if(contains(prevNode) && !effectiveNodeSet.contains(prevNode))
						{
							effectiveNodeSet.add(prevNode);
							newFreshNodes.add(prevNode);
						}
					}				
				}
			}
			curFreshNodes = newFreshNodes;
		}
		
		ArrayList<T> intermidiateNodeList = new ArrayList<T>(effectiveNodeSet);
		Collections.sort(intermidiateNodeList);
		return intermidiateNodeList;
	}

	protected void makeAllPossibleMasses(boolean makeEdgeMap)
	{
		HashSet<T> nodes = new HashSet<T>();
		ArrayList<T> curFreshNodes = new ArrayList<T>();
		
		T zero = getZero();
		nodes.add(zero);
		curFreshNodes.add(zero);

		if(makeEdgeMap)
		{
			edgeMap = new HashMap<T, ArrayList<DeNovoGraph.Edge<T>>>();
			edgeMap.put(zero, new ArrayList<DeNovoGraph.Edge<T>>());
		}
		
		boolean isSource = true;
		for(int i=0; i<maxLength; i++)
		{
			ArrayList<T> newFreshNodes = new ArrayList<T>();
			for(T node : curFreshNodes)
			{
				for(AminoAcid aa : aaSet)
				{
					T newNode = getNextNode(node, aa);
					assert(newNode != null): node.getNominalMass()+" "+aa.getResidueStr();
					DeNovoGraph.Edge<T> edge = new DeNovoGraph.Edge<T>(node, aa.getProbability(), aaSet.getIndex(aa), aa.getMass());
					if(isSource && enzyme != null)
					{
						if(enzyme.isCleavable(aa))
							edge.setCleavageScore(enzyme.getPeptideCleavageCredit());
						else
							edge.setCleavageScore(enzyme.getPeptideCleavagePenalty());
					}
					if(nodes.add(newNode))	// newly generated node
					{
						ArrayList<DeNovoGraph.Edge<T>> edges = new ArrayList<DeNovoGraph.Edge<T>>();
						edges.add(edge);
						if(makeEdgeMap)
							edgeMap.put(newNode, edges);
						
						newFreshNodes.add(newNode);
					}
					else	// existing node
					{
						if(makeEdgeMap)
							edgeMap.get(newNode).add(edge);
					}
				}
				if(isSource)
					isSource = false;
			}
			curFreshNodes = newFreshNodes;
		}
		
		allNodes = new ArrayList<T>(nodes);
		Collections.sort(allNodes);
	}

	@Override
	public Sequence<T> toCumulativeSequence(boolean isPrefix, Peptide pep) 
	{
		Sequence<T> cumSeq = new Sequence<T>();
		
		T curNode = getZero();
		for(int i=pep.size()-1; i>=0; i--)
		{
			AminoAcid aa;
			if(isPrefix)
				aa = pep.get(pep.size()-1-i);
			else
				aa = pep.get(i);
			curNode = getNextNode(curNode, aa);
			cumSeq.add(curNode);
		}
		return cumSeq;
	}
}
