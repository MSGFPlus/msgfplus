package msgf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import suffixarray.SuffixArray;

import msutil.*;


public class GeneratingFunction<T extends Matter> implements GF<T> {
	private final DeNovoGraph<T> graph;
	
	private boolean backtrack = true;
	private boolean calcNumber = true;
	private boolean calcProb = true;
	private Enzyme enzyme = Enzyme.TRYPSIN;
	
	private int numScoreBinsPerNode = 1000;
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
	public GeneratingFunction<T> numScoreBinsPerNode(int numBins)	{ this.numScoreBinsPerNode = numBins; return this; }
	public GeneratingFunction<T> gfTableCapacity(int gfTableCapacity) { this.gfTableCapacity = gfTableCapacity; return this;}
	
	public boolean backtrack()		{ return backtrack; }
	public boolean calcNumber()		{ return calcNumber; }
	public boolean calcProb()		{ return calcProb; }
	public Enzyme getEnzyme()	{ return enzyme; }
	public int getNumScoreBinsPerNode()	{ return numScoreBinsPerNode; }
	public boolean isGFComputed()	{ return this.isGFComputed; }
	public DeNovoGraph<T> getGraph()	{ return graph; }
	
	protected HashMap<T, ScoreDist> getFwdTable()	{ return fwdTable; }
	protected BacktrackTable<T> getBacktrackTable()	{ return backtrackTable; }
//	public HashMap<T, Integer> getSRMScore()	{ return srmScore; }
//	//added by kyowon
//	public void setSRMScore(HashMap<T, Integer> srmScore) { this.srmScore = srmScore; }
	
	@Override
	public int getScore(Annotation annotation)
	{
		return graph.getScore(annotation);
	}
	
	public int getEnergy(Annotation annotation) {
		return getMaxScore()-getScore(annotation);
	}
	
	public float getSpectralProbability(Annotation annotation) {
		int score = getScore(annotation);
		return getSpectralProbability(score);
	}
	
	// score: inclusive
	@Override
	public float getSpectralProbability(int score) {
		if(!this.distribution.isProbSet())
			return 100;
		return distribution.getSpectralProbability(score);
	}
	
	public float getNumEqualBetterPeptides(Annotation annotation)
	{
		int score = getScore(annotation);
		return getNumEqualOrBetterPeptides(score);
	}

	public float getNumEqualOrBetterPeptides(int score)
	{
		if(!this.distribution.isNumSet())
			return -1;
		return distribution.getNumEqualOrBetterPeptides(score);
	}

	public float getDictionarySize(float specProb)
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
	
	@Override
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
				backtrackTable.getReconstructions(sink, score-enzyme.getNeighboringAACleavageCredit(), "R.", reconstructions, sa);
			for(T sink : graph.getSinkList())
				backtrackTable.getReconstructions(sink, score-enzyme.getNeighboringAACleavagePenalty(), "L.", reconstructions, sa);
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
	
	@Override
	public int getMaxScore()	{ return this.distribution.getMaxScore(); }
	
	@Override
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
		if(enzyme != null)
		{
			int neighboringAACleavageCredit = enzyme.getNeighboringAACleavageCredit();
			int neighboringAACleavagePenalty = enzyme.getNeighboringAACleavagePenalty();
			finalDist = factory.getInstance(mergedDist.getMinScore()+neighboringAACleavagePenalty, mergedDist.getMaxScore()+neighboringAACleavageCredit);
			if(calcNumber)
			{
				finalDist.addNumDist(mergedDist, neighboringAACleavageCredit, enzyme.getResidues().size());
				finalDist.addNumDist(mergedDist, neighboringAACleavagePenalty, graph.getAASet().size()-enzyme.getResidues().size());
			}
			if(calcProb)
			{
				finalDist.addProbDist(mergedDist, neighboringAACleavageCredit, enzyme.getProbCleavageSites());
				finalDist.addProbDist(mergedDist, neighboringAACleavagePenalty, 1-enzyme.getProbCleavageSites());
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
		int curMaxScore = Integer.MIN_VALUE;;
		int curMinScore = Integer.MAX_VALUE;
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
				if(prevDist.getMinScore()+combinedScore < curMinScore)
					curMinScore = prevDist.getMinScore() + combinedScore;
				edges.add(edge);
			}
		}
		if(curMinScore >= curMaxScore)
			return;
		if(curMaxScore-curMinScore > numScoreBinsPerNode)
			curMinScore = curMaxScore-numScoreBinsPerNode;

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
				curDist.setProb(curDist.maxScore-1, Float.MIN_VALUE);
		}
		fwdTable.put(curNode, curDist);
		if(backtrack)
			backtrackTable.put(curNode, backPointer);	
	}

}
