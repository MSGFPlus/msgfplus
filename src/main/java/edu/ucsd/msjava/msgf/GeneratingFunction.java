package edu.ucsd.msjava.msgf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;


import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.suffixarray.SuffixArray;


public class GeneratingFunction<T extends Matter> implements GF<T> {
	private final DeNovoGraph<T> graph;
	
	private boolean backtrack = true;
	private boolean calcNumber = true;
	private boolean calcProb = true;
	private Enzyme enzyme = Enzyme.TRYPSIN;
	
//	private int numScoreBinsPerNode = 1000;
	private int gfTableCapacity;
	
	private ScoreDist distribution = null;
	private BacktrackTable<T> backtrackTable = null;	
	
	private class GFTable extends LinkedHashMap<T,ScoreDist> {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private final int capacity;
		public GFTable(int capacity)
		{
			super(capacity+1, 1.1f, false);
			this.capacity = capacity;
		}
		
		@Override
		protected boolean removeEldestEntry(Map.Entry<T,ScoreDist> eldest)
		{
			return size() > capacity;
		}
	}
	
	private HashMap<T, ScoreDist> fwdTable;
	private HashMap<T, Integer> minScoreTable = null;
	
	private boolean isGFComputed = false;
	
//	private HashMap<T, Integer> srmScore = null;
	
	public GeneratingFunction(DeNovoGraph<T> graph)
	{
		this.graph = graph;
		this.gfTableCapacity = 1+graph.intermediateNodes.size()+graph.sinkNodes.size();
	}
	
	// Builder
	public GeneratingFunction<T> doNotBacktrack()	{ this.backtrack = false; return this; }
	public GeneratingFunction<T> doNotCalcNumber()	{ this.calcNumber = false; return this; }
	public GeneratingFunction<T> doNotCalcProb()	{ this.calcProb = false; return this; }
	public GeneratingFunction<T> enzyme(Enzyme enzyme)	{ this.enzyme = enzyme; return this; }
//	public GeneratingFunction<T> numScoreBinsPerNode(int numBins)	{ this.numScoreBinsPerNode = numBins; return this; }
	public GeneratingFunction<T> gfTableCapacity(int gfTableCapacity) { this.gfTableCapacity = gfTableCapacity; return this;}
	
	public boolean backtrack()		{ return backtrack; }
	public boolean calcNumber()		{ return calcNumber; }
	public boolean calcProb()		{ return calcProb; }
	public Enzyme getEnzyme()	{ return enzyme; }
//	public int getNumScoreBinsPerNode()	{ return numScoreBinsPerNode; }
	public boolean isGFComputed()	{ return this.isGFComputed; }
	public DeNovoGraph<T> getGraph()	{ return graph; }
	
	protected HashMap<T, ScoreDist> getFwdTable()	{ return fwdTable; }
	protected BacktrackTable<T> getBacktrackTable()	{ return backtrackTable; }
	
	public int getScore(Annotation annotation)
	{
		return graph.getScore(annotation);
	}
	
	public int getEnergy(Annotation annotation) {
		return getMaxScore()-getScore(annotation);
	}
	
	public double getSpectralProbability(Annotation annotation) {
		int score = getScore(annotation);
		return getSpectralProbability(score);
	}
	
	// score: inclusive
	public double getSpectralProbability(int score) {
		if(!this.distribution.isProbSet())
			return 100;
		return distribution.getSpectralProbability(score);
	}
	
	public double getNumEqualBetterPeptides(Annotation annotation)
	{
		int score = getScore(annotation);
		return getNumEqualOrBetterPeptides(score);
	}

	public double getNumEqualOrBetterPeptides(int score)
	{
		if(!this.distribution.isNumSet())
			return -1.;
		return distribution.getNumEqualOrBetterPeptides(score);
	}

	public double getDictionarySize(float specProb)
	{
		return getNumEqualOrBetterPeptides(getThresholdScore(specProb));
	}
	
	// returns t where totalProb(t) > specProb && totalProb(t+1) <= specProb
	public static int getThresholdScore(float specProb, ScoreDist distribution)
	{
		if(!distribution.isProbSet())
			return -1;
		float totalProb = 0;
		
		for(int t = distribution.getMaxScore()-1; t>=distribution.getMinScore(); t--)
		{
			totalProb += distribution.getProbability(t);
			if(totalProb > specProb)
				return t;
		}
		return -1;
	}
	
	// returns t where totalProb(t) > specProb && totalProb(t+1) <= specProb
	public int getThresholdScore(float specProb)
	{
		return getThresholdScore(specProb, distribution);
	}
	
	public ScoreDist getScoreDist()
	{
		return distribution;
	}
	
	/**
	 * Generate reconstructions with score "score" and have match with "sa" and put it in "reconstructions".
	 * @param score the score of reconstructions to be generated
	 * @param reconstructions a container where reconstructions will be stored
	 * @param sa suffix array that will filter reconstructions
	 * @return
	 */
	private void generateReconstructions(int score, ArrayList<String> reconstructions, SuffixArray sa)
	{
		if(backtrackTable == null)
			return;
		if(enzyme == null)
		{
			for(T sink : graph.getSinkList())
				backtrackTable.getReconstructions(sink, score, "", reconstructions, sa);
		}
		else
		{
			//TODO: add prefix info?
			for(T sink : graph.getSinkList())
				backtrackTable.getReconstructions(sink, score-graph.getAASet().getNeighboringAACleavageCredit(), "R.", reconstructions, sa);
			for(T sink : graph.getSinkList())
				backtrackTable.getReconstructions(sink, score-graph.getAASet().getNeighboringAACleavagePenalty(), "L.", reconstructions, sa);
		}
	}

	public String getOneReconstruction(int score)
	{
		if(backtrackTable == null)
			return null;
		return backtrackTable.getOneReconstruction(graph.getPMNode(), score, "");
	}
	
	public ArrayList<String> getReconstructions(int score)
	{
		ArrayList<String> reconstructions = new ArrayList<String>();
		generateReconstructions(score, reconstructions, null);
		return reconstructions;
	}
	
	public ArrayList<String> getReconstructionsEqualOrAboveScore(int score)
	{
		ArrayList<String> reconstructions = new ArrayList<String>();
		for(int t=this.getMaxScore()-1; t>=score; t--)
			generateReconstructions(t, reconstructions, null);
		return reconstructions;
	}
	
	public ArrayList<String> getDictionary(float specProbThreshold)
	{
		assert(calcProb);
		int threshold = getThresholdScore(specProbThreshold);
		return getReconstructionsEqualOrAboveScore(threshold+1);
	}

	public ArrayList<String> getReconstructions(float specProbThreshold, float numRecsThreshold, boolean isNumInclusive, SuffixArray sa)
	{
		assert(calcProb && calcNumber);
		ArrayList<String> recs = new ArrayList<String>();
		int threshold = getThresholdScore(specProbThreshold);
		float numRecs = 0;
		for(int t = getMaxScore()-1; t>threshold; t--)
		{
			numRecs += distribution.getNumberRecs(t);
			if(!isNumInclusive)
			{
				if(numRecs <= numRecsThreshold)
					generateReconstructions(t, recs, sa);
				else
					break;
			}
			else
			{
				generateReconstructions(t, recs, sa);
				if(numRecs >= numRecsThreshold)
					break;
			}
		}
		return recs;
	}
	
	public int getMinScore()	{ return this.distribution.getMinScore(); }
	
	public int getMaxScore()	{ return this.distribution.getMaxScore(); }
	
	public void setUpScoreThreshold(int score)
	{
		minScoreTable = new HashMap<T, Integer>();
		if(enzyme != null)
			score -= graph.getAASet().getNeighboringAACleavageCredit();
		
		for(T sink : graph.getSinkList())
		{
			minScoreTable.put(sink, score);
			for(DeNovoGraph.Edge<T> edge : graph.getEdges(sink))
			{
				T prevNode = edge.getPrevNode();
				int newPrevMinScore = score - edge.getEdgeScore();
				Integer prevMinScore = minScoreTable.get(prevNode);
				if(prevMinScore == null || prevMinScore > newPrevMinScore)
					minScoreTable.put(prevNode, newPrevMinScore);
			}
		}
		
		ArrayList<T> intermediateNodeList = graph.getIntermediateNodeList();
		
		for(int i=intermediateNodeList.size()-1; i>=0; i--)
		{
			T curNode = intermediateNodeList.get(i);
			Integer curScore = minScoreTable.get(curNode);
			if(curScore == null)
				continue;
			int curNodeScore = graph.getNodeScore(curNode);
			for(DeNovoGraph.Edge<T> edge : graph.getEdges(curNode))
			{
				T prevNode = edge.getPrevNode();
				int newPrevMinScore = curScore - (curNodeScore + edge.getEdgeScore());
				Integer prevMinScore = minScoreTable.get(prevNode);
				if(prevMinScore == null || prevMinScore > newPrevMinScore)
					minScoreTable.put(prevNode, newPrevMinScore);
			}
		}
	}
	
	public boolean computeGeneratingFunction()
	{
		ScoreDistFactory factory = new ScoreDistFactory(calcNumber, calcProb);
		// initialization of the source
		ScoreDist sourceDist = factory.getInstance(0, 1);
		if(calcNumber)
			sourceDist.setNumber(0, 1);
		if(calcProb)	
			sourceDist.setProb(0, 1);
		fwdTable = new GFTable(gfTableCapacity);
		fwdTable.put(graph.getSource(), sourceDist);
		if(backtrack)
		{
			backtrackTable = new BacktrackTable<T>(graph);
			BacktrackPointer sourcePointer = new BacktrackPointer(0, 1, 0);
			sourcePointer.setBacktrack(0, 0);
			backtrackTable.put(graph.getSource(), sourcePointer);
		}
		
		// dynamic programming, source node (i=0) is excluded
		ArrayList<T> intermediateNodeList = graph.getIntermediateNodeList();
		
		for(int i=1; i<intermediateNodeList.size(); i++)
		{
			T curNode = intermediateNodeList.get(i);
			setCurNode(curNode, factory);
		}
		
		// process dest node
		int minScore = Integer.MAX_VALUE;
		int maxScore = Integer.MIN_VALUE;
		
		for(T curNode : graph.getSinkList())
		{
			setCurNode(curNode, factory);
			ScoreDist curDist = fwdTable.get(curNode);
			if(curDist == null)	// curNode is not connected from the source
				continue;
			if(curDist.getMinScore() < minScore)
				minScore = curDist.getMinScore();
			if(curDist.getMaxScore() > maxScore)
				maxScore = curDist.getMaxScore();
		}		

		if(maxScore <= minScore)
			return false;
		
		// merge distributions of dest nodes
		ScoreDist mergedDist = factory.getInstance(minScore, maxScore);
		for(T sinkNode : graph.getSinkList())
		{
			if(calcNumber)
				mergedDist.addNumDist(fwdTable.get(sinkNode), 0);
			if(calcProb)
				mergedDist.addProbDist(fwdTable.get(sinkNode), 0, 1);
		}
		
		// process neighboring amino acid
		ScoreDist finalDist;
		if(enzyme != null && enzyme.getResidues() != null)
		{
			int neighboringAACleavageCredit = graph.getAASet().getNeighboringAACleavageCredit();
			int neighboringAACleavagePenalty = graph.getAASet().getNeighboringAACleavagePenalty();
			finalDist = factory.getInstance(mergedDist.getMinScore()+neighboringAACleavagePenalty, mergedDist.getMaxScore()+neighboringAACleavageCredit);
			if(calcNumber)
			{
				finalDist.addNumDist(mergedDist, neighboringAACleavageCredit, enzyme.getResidues().length);
				finalDist.addNumDist(mergedDist, neighboringAACleavagePenalty, graph.getAASet().size()-enzyme.getResidues().length);
			}
			if(calcProb)
			{
				finalDist.addProbDist(mergedDist, neighboringAACleavageCredit, graph.getAASet().getProbCleavageSites());
				finalDist.addProbDist(mergedDist, neighboringAACleavagePenalty, 1-graph.getAASet().getProbCleavageSites());
			}
		}
		else
		{
			finalDist = mergedDist;
		}
		
		this.distribution = finalDist;
		isGFComputed = true;
		return true;
	}

	// scoreThreshold : inclusive
	public HashMap<T, Float> getDestProfile(int scoreThreshold)
	{
		assert(calcNumber);
		HashMap<T, Float> destProf = new HashMap<T, Float>();
		for(T sinkNode : graph.getSinkList())
		{
			float num = 0;
			ScoreDist dist = fwdTable.get(sinkNode);
			for(int t = dist.getMaxScore()-1; t >= dist.getMinScore() && t>=scoreThreshold; t--)
				num += dist.getNumberRecs(t);
			if(num > 0)
				destProf.put(sinkNode, num);
		}
		return destProf;
	}
	
	private void setCurNode(T curNode, ScoreDistFactory scoreDistFactory)
	{
		int curNodeScore = graph.getNodeScore(curNode);
		int curMaxScore = Integer.MIN_VALUE;
		int curMinScore;
		if(minScoreTable == null)
			curMinScore = Integer.MAX_VALUE;
		else
		{
			Integer min = minScoreTable.get(curNode);
			if(min == null)
				return;
			curMinScore = min;
		}
		// determine minScore and maxScore
		ArrayList<DeNovoGraph.Edge<T>> edges = new ArrayList<DeNovoGraph.Edge<T>>(); // modified by kyowon
		for(DeNovoGraph.Edge<T> edge : graph.getEdges(curNode))
		{
			T prevNode = edge.getPrevNode();
			ScoreDist prevDist = fwdTable.get(prevNode);
			if(prevDist != null)
			{
				int edgeScore = edge.getEdgeScore();
				int combinedScore = curNodeScore + edgeScore;
				if(prevDist.getMaxScore()+combinedScore > curMaxScore)
					curMaxScore = prevDist.getMaxScore() + combinedScore;
				if(minScoreTable == null)
				{
					if(prevDist.getMinScore()+combinedScore < curMinScore)
						curMinScore = prevDist.getMinScore() + combinedScore;
				}
				edges.add(edge);
			}
		}
		if(curMinScore >= curMaxScore)
			return;

		ScoreDist curDist = scoreDistFactory.getInstance(curMinScore, curMaxScore);
		BacktrackPointer backPointer = null;
		if(backtrack)
			backPointer = new BacktrackPointer(curMinScore, curMaxScore, curNodeScore);
		for(DeNovoGraph.Edge<T> edge : edges)
		{
			T prevNode = edge.getPrevNode();
			ScoreDist prevDist = fwdTable.get(prevNode);
			if(prevDist != null)
			{
				int edgeScore = edge.getEdgeScore();
				int combinedScore = curNodeScore + edgeScore;

				if(calcNumber)
					curDist.addNumDist(prevDist, combinedScore, 1);
				if(calcProb)
					curDist.addProbDist(prevDist, combinedScore, edge.getEdgeProbability());
				if(backtrack)
				{
					BacktrackPointer prevPointer = backtrackTable.get(prevNode);
					backPointer.addBacktrackPointers(prevPointer, edge.getEdgeIndex(), edgeScore);
				}
			}
		}
		if(calcProb)
		{
			if(curDist.getProbability(curDist.maxScore-1) == 0)	// to avoid underflow
			{
				assert(false): "Underflow! " + curNode.getNominalMass()+ " " + curDist.getProbability(curDist.maxScore-1);
				curDist.setProb(curDist.maxScore-1, Float.MIN_VALUE);
			}
		}
		fwdTable.put(curNode, curDist);
		if(backtrack)
			backtrackTable.put(curNode, backPointer);	
	}
}
