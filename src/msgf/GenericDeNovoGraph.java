package msgf;

import java.util.ArrayList;
import java.util.HashMap;

import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Annotation;
import msutil.Composition;
import msutil.Enzyme;
import msutil.Matter;
import msutil.Modification;
import msutil.ModifiedAminoAcid;
import msutil.Peptide;
import msutil.Sequence;

public class GenericDeNovoGraph<T extends Matter> extends DeNovoGraph<T> {
	private HashMap<T, Integer> nodeScore;
	private ScoredSpectrum<T> scoredSpec;
	private DeNovoNodeFactory<T> factory;
	private Enzyme enzyme;
	private HashMap<T, ArrayList<DeNovoGraph.Edge<T>>> edgeMap;
	
	public GenericDeNovoGraph(DeNovoNodeFactory<T> factory, float parentMass, Tolerance pmTolerance, Enzyme enzyme, ScoredSpectrum<T> scoredSpec) 
	{
		this.factory = factory;
		this.enzyme = enzyme;
		this.scoredSpec = scoredSpec;
		edgeMap = new HashMap<T, ArrayList<DeNovoGraph.Edge<T>>>();

		float peptideMass = parentMass - (float)Composition.H2O;
		super.sinkNodes = factory.getNodes(peptideMass, pmTolerance);
		setEdgesToSinkNodes();
		
		super.intermediateNodes = factory.getIntermediateNodes(super.sinkNodes);
		super.pmNode = null;
		if(sinkNodes.size() > 0)
		{
			super.source = super.intermediateNodes.get(0);
			super.pmNode = sinkNodes.get(0);
			for(int i=1; i<sinkNodes.size(); i++)
				if(Math.abs(peptideMass - sinkNodes.get(i).getMass()) < Math.abs(peptideMass - pmNode.getMass()))
					pmNode = sinkNodes.get(i);
		}
		super.destination = factory.getInfinity();
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
					edgeScore += enzyme.getPeptideCleavageCredit();
				else
					edgeScore += enzyme.getPeptideCleavagePenalty();
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
				score += enzyme.getNeighboringAACleavageCredit();
			else
				score += enzyme.getNeighboringAACleavagePenalty();
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

	private void setEdgesToSinkNodes()
	{
		boolean isReverse = factory.isReverse();
		for(T curNode : sinkNodes)
		{
			ArrayList<DeNovoGraph.Edge<T>> unmodEdges = factory.getEdges(curNode);
			AminoAcidSet aaSet = factory.getAASet();
			
			// apply fixed mods
			
			// variable mods
			ArrayList<DeNovoGraph.Edge<T>> modEdges = new ArrayList<DeNovoGraph.Edge<T>>(); 
			ArrayList<ModifiedAminoAcid> termVariableMods;
			if(isReverse)	// consider N-term specific mods
				termVariableMods = aaSet.getNTermVariableMods();
			else			// consider c-term specific mods
				termVariableMods = aaSet.getCTermVariableMods();
			
			// non residue specific mods
			for(ModifiedAminoAcid modAA : termVariableMods)
			{
				
			}
			
			for(ModifiedAminoAcid modAA : termVariableMods)
			{
				if(modAA.isResidueSpecific())
				{
					T prevNode = factory.getPreviousNode(curNode, modAA);
					if(prevNode != null)
					{
						DeNovoGraph.Edge<T> edge = new DeNovoGraph.Edge<T>(prevNode, modAA.getProbability(), aaSet.getIndex(modAA), modAA.getMass());
						modEdges.add(edge);
					}
				}
			}
			for(DeNovoGraph.Edge<T> edge : unmodEdges)
				edge.setErrorScore(0);
			edgeMap.put(curNode, unmodEdges);
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
}
