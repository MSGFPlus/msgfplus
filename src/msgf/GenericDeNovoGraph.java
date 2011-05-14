package msgf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Annotation;
import msutil.Composition;
import msutil.Enzyme;
import msutil.Matter;
import msutil.Peptide;
import msutil.Sequence;
import msutil.Modification.Location;

public class GenericDeNovoGraph<T extends Matter> extends DeNovoGraph<T> {
	private HashMap<T, Integer> nodeScore;
	private ScoredSpectrum<T> scoredSpec;
	private DeNovoNodeFactory<T> factory;
	private Enzyme enzyme;
	private HashMap<T, ArrayList<DeNovoGraph.Edge<T>>> edgeMap;
	
	public GenericDeNovoGraph(
			DeNovoNodeFactory<T> factory, 
			float parentMass, 
			Tolerance pmTolerance, 
			Enzyme enzyme, 
			ScoredSpectrum<T> scoredSpec
		) 
	{
		this(factory, parentMass, pmTolerance, enzyme, scoredSpec, false);
	}
	
	public GenericDeNovoGraph(
			DeNovoNodeFactory<T> factory, 
			float parentMass, 
			Tolerance pmTolerance, 
			Enzyme enzyme, 
			ScoredSpectrum<T> scoredSpec,
			boolean containsModifiedSinkEdge
		) 
	{
		this.factory = factory;
		this.enzyme = enzyme;
		this.scoredSpec = scoredSpec;
		edgeMap = new HashMap<T, ArrayList<DeNovoGraph.Edge<T>>>();

		super.source = factory.getZero();
		
		float peptideMass = parentMass - (float)Composition.H2O;
		super.pmNode = factory.getNode(peptideMass);
		
		super.sinkNodes = factory.getNodes(peptideMass, pmTolerance);
		setEdgesToSinkNodes(containsModifiedSinkEdge);
		
		setIntermediateNodes();
		setEdgesToIntermediateNodes();
		computeNodeScores();
	}

	@Override
	public T getComplementNode(T srm) {
		return factory.getComplementNode(srm, pmNode);
	}
	
	@Override
	public ArrayList<DeNovoGraph.Edge<T>> getEdges(T curNode) 
	{
		return edgeMap.get(curNode);
	}

	@Override
	public int getNodeScore(T node) {
		return nodeScore.get(node);
	}

	@Override
	public int getScore(Peptide pep) {
		int score = 0;
		
		T prevNode = source;
//		score += getNodeScore(source);
		Sequence<T> seq = factory.toCumulativeSequence(!isReverse(), pep);
		T pmNode = seq.get(seq.size()-1);
		if(!sinkNodes.contains(pmNode))
			return Integer.MIN_VALUE;
		for(int i=0; i<seq.size()-1; i++)
		{
			T curNode = seq.get(i);
			int nodeScore = getNodeScore(curNode);
			
			AminoAcid aa;
			if(!isReverse())
				aa = pep.get(i);
			else
				aa = pep.get(pep.size()-1-i);
			/////////// must be implemented again
			int edgeScore = scoredSpec.getEdgeScore(curNode, prevNode, aa.getMass());
			if(prevNode == source && enzyme != null)
			{
				if(enzyme.isCleavable(aa))
					edgeScore += factory.getAASet().getPeptideCleavageCredit();
				else
					edgeScore += factory.getAASet().getPeptideCleavagePenalty();
			}
			////////////////
			score += nodeScore + edgeScore;
//			System.out.println(curNode.getNominalMass()+" "+nodeScore+" "+edgeScore+" "+score);
			prevNode = curNode;
		}
		return score;
	}

	@Override
	public int getScore(Annotation annotation) {
		int score = getScore(annotation.getPeptide());
		if(score > Integer.MIN_VALUE && enzyme != null)
		{
			AminoAcid neighboringAA;
			if(enzyme.isCTerm())
				neighboringAA = annotation.getPrevAA();
			else
				neighboringAA = annotation.getNextAA();
			if(neighboringAA == null || enzyme.isCleavable(neighboringAA))
				score += factory.getAASet().getNeighboringAACleavageCredit();
			else
				score += factory.getAASet().getNeighboringAACleavagePenalty();
		}
		
		return score;
	}

	@Override
	public boolean isReverse() {
		return factory.isReverse();
	}

	@Override
	public AminoAcidSet getAASet() {
		return factory.getAASet();
	}
	
	private void computeNodeScores() {
		nodeScore = new HashMap<T,Integer>();
		boolean isReverse = factory.isReverse();
		
		nodeScore.put(source, 0);
		for(int i=1; i<intermediateNodes.size(); i++)
		{
			T node = intermediateNodes.get(i);
			T compNode = this.getComplementNode(node);
			int score;
			if(isReverse)
				score = scoredSpec.getNodeScore(compNode, node);
			else
				score = scoredSpec.getNodeScore(node, compNode);
			nodeScore.put(node, score);
		}
		for(T node : this.sinkNodes)
			nodeScore.put(node, 0);
	}

	private void setEdgesToSinkNodes(boolean containsModifiedSinkEdge)
	{
		boolean isReverse = factory.isReverse();
		Location location;
		if(isReverse)
		{
			if(!containsModifiedSinkEdge)
				location = Location.N_Term;
			else
				location = Location.Protein_N_Term;
		}
		else
		{
			if(!containsModifiedSinkEdge)
				location = Location.C_Term;
			else
				location = Location.Protein_C_Term;
		}
		
		AminoAcidSet aaSet = factory.getAASet();
		ArrayList<AminoAcid> aaList = aaSet.getAAList(location);
		for(T curNode : sinkNodes)
		{
			ArrayList<DeNovoGraph.Edge<T>> edges = new ArrayList<DeNovoGraph.Edge<T>>();
			for(AminoAcid aa : aaList)
			{
				T prevNode = factory.getPreviousNode(curNode, aa);
				if(prevNode != null)
				{
					DeNovoGraph.Edge<T> edge = new DeNovoGraph.Edge<T>(prevNode, aa.getProbability(), aaSet.getIndex(aa), aa.getMass());				
					edges.add(edge);
				}
			}
			edgeMap.put(curNode, edges);
		}		
	}
	
	private void setEdgesToIntermediateNodes()
	{
		for(T curNode : intermediateNodes)
		{
			ArrayList<DeNovoGraph.Edge<T>> edges = factory.getEdges(curNode);
			for(DeNovoGraph.Edge<T> edge : edges)
			{
				T prevNode = edge.getPrevNode();
				int errorScore = scoredSpec.getEdgeScore(curNode, prevNode, edge.getEdgeMass());
				edge.setErrorScore(errorScore);
			}
			edgeMap.put(curNode, edges);
		}
	}	
	
	private void setIntermediateNodes() 
	{
		HashSet<T> depth1Nodes = new HashSet<T>();
		for(T sink : sinkNodes)
		{
			ArrayList<DeNovoGraph.Edge<T>> edges = getEdges(sink);
			if(edges != null)
			{
				for(DeNovoGraph.Edge<T> edge : edges)
				{
					T prevNode = edge.getPrevNode();
					depth1Nodes.add(prevNode);
				}				
			}
		}
		
		ArrayList<T> intermidiateNodeList = factory.getLinkedNodeList(depth1Nodes);
		Collections.sort(intermidiateNodeList);
		this.intermediateNodes = intermidiateNodeList;
	}
	
}
