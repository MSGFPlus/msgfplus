package edu.ucsd.msjava.msgf;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.Matter;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.Sequence;
import edu.ucsd.msjava.msutil.Modification.Location;

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
	public boolean isReverse()
	{
		return enzyme == null || enzyme.isCTerm();
	}
	
	public int getMaxLength()	{ return maxLength; }
	
	public ArrayList<T> getAllNodes()	{ return allNodes; }
	
	public int size() {
		return allNodes.size(); 
	}
	
	public AminoAcidSet getAASet() {
		return aaSet;
	}

	public DeNovoGraph.Edge<T> getEdge(T curNode, T prevNode)
	{
		for(DeNovoGraph.Edge<T> edge : getEdges(curNode))
		{
			if(edge.getPrevNode().equals(prevNode))
				return edge;
		}
		return null;
	}
	
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
	
	public abstract T getZero();
	
	public Enzyme getEnzyme()
	{
		return enzyme;
	}
	
	public ArrayList<T> getLinkedNodeList(Collection<T> destNodes) 
	{
		HashSet<T> effectiveNodeSet = new HashSet<T>(destNodes);
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
		
		T zero = getZero();
		nodes.add(zero);

		if(makeEdgeMap)
		{
			edgeMap = new HashMap<T, ArrayList<DeNovoGraph.Edge<T>>>();
			edgeMap.put(zero, new ArrayList<DeNovoGraph.Edge<T>>());
		}
		
		// length 1
		ArrayList<T> curFreshNodes = new ArrayList<T>();
		Location location;
		if(isReverse())	// C-term
			location = Location.C_Term;
		else			// N-term
			location = Location.N_Term;
		
		for(AminoAcid aa : aaSet.getAAList(location))
		{
			T newNode = getNextNode(zero, aa);
			boolean isNewNode = nodes.add(newNode);
			if(isNewNode)
				curFreshNodes.add(newNode);
			
			if(makeEdgeMap)
			{
				DeNovoGraph.Edge<T> edge = new DeNovoGraph.Edge<T>(zero, aa.getProbability(), aaSet.getIndex(aa), aa.getMass());
				if(enzyme != null)
				{
					if(enzyme.isCleavable(aa))
						edge.setCleavageScore(aaSet.getPeptideCleavageCredit());
					else
						edge.setCleavageScore(aaSet.getPeptideCleavagePenalty());
				}
				if(isNewNode)	// newly generated node
				{
					ArrayList<DeNovoGraph.Edge<T>> edges = new ArrayList<DeNovoGraph.Edge<T>>();
					edges.add(edge);
					edgeMap.put(newNode, edges);
				}
				else	// existing node
				{
					edgeMap.get(newNode).add(edge);
				}
			}
		}
		
		// length >=2
		for(int i=1; i<maxLength; i++)
		{
			ArrayList<T> newFreshNodes = new ArrayList<T>();
			for(T node : curFreshNodes)
			{
				for(AminoAcid aa : aaSet)
				{
					T newNode = getNextNode(node, aa);
					assert(newNode != null): node.getNominalMass()+" "+aa.getResidueStr();
					boolean isNewNode = nodes.add(newNode);
					if(isNewNode)
						newFreshNodes.add(newNode);
					if(makeEdgeMap)
					{
						DeNovoGraph.Edge<T> edge = new DeNovoGraph.Edge<T>(node, aa.getProbability(), aaSet.getIndex(aa), aa.getMass());
						if(isNewNode)	// newly generated node
						{
							ArrayList<DeNovoGraph.Edge<T>> edges = new ArrayList<DeNovoGraph.Edge<T>>();
							edges.add(edge);
							edgeMap.put(newNode, edges);
						}
						else	// existing node
						{
							edgeMap.get(newNode).add(edge);
						}
					}
				}
			}
			curFreshNodes = newFreshNodes;
		}
		
		allNodes = new ArrayList<T>(nodes);
		Collections.sort(allNodes);
	}

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
