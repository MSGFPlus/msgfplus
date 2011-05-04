package msgf2d;

import java.util.ArrayList;
import java.util.Hashtable;

import suffixarray.SuffixArray;

import msgf.DeNovoGraph;
import msgf.ScoredSpectrum;
import msgf.DeNovoGraph.Edge;
import msutil.AminoAcidSet;
import msutil.Annotation;
import msutil.Matter;
import msutil.Peptide;

public class GeneratingFunction2D<T extends Matter> {
	private final ScoredSpectrum<T> scoredSpec1;
	private final ScoredSpectrum<T> scoredSpec2;
	private final DeNovoGraph<T> graph;
	private final AminoAcidSet aaSet;
	
	public GeneratingFunction2D(ScoredSpectrum<T> scoredSpec1, ScoredSpectrum<T> scoredSpec2, DeNovoGraph<T> graph)
	{
		this.scoredSpec1 = scoredSpec1;
		this.scoredSpec2 = scoredSpec2;
		this.graph = graph;
		this.aaSet = graph.getAASet();
	}
	
	private boolean backtrack = true;
	private boolean calcNumber = true;
	private boolean calcProb = true;
	
	private int numScoreBinsPerNode = 1000;
	
	private ScoreDistMerged distribution = null;
	private BacktrackTable2D<T> backtrackTable = null;	
	private Hashtable<T, ScoreDist2D> fwdTable;
	
	private Hashtable<T, Integer> srmScore1 = null;	
	private Hashtable<T, Integer> srmScore2 = null;	
	
	// Builder
	public GeneratingFunction2D<T> doNotBacktrack()	{ this.backtrack = false; return this; }
	public GeneratingFunction2D<T> doNotCalcNumber()	{ this.calcNumber = false; return this; }
	public GeneratingFunction2D<T> doNotCalcProb()	{ this.calcProb = false; return this; }
	public GeneratingFunction2D<T> numScoreBinsPerNode(int numBins)	{ this.numScoreBinsPerNode = numBins; return this; }
	
	public boolean backtrack()		{ return backtrack; }
	public boolean calcNumber()		{ return calcNumber; }
	public boolean calcProb()		{ return calcProb; }
	public int getNumScoreBinsPerNode()	{ return numScoreBinsPerNode; }
	
	protected Hashtable<T, ScoreDist2D> getFwdTable()	{ return fwdTable; }
	protected DeNovoGraph<T> getGraph()	{ return graph; }
	protected ScoredSpectrum<T> getScoredSpectrum1()	{ return scoredSpec1; }
	protected ScoredSpectrum<T> getScoredSpectrum2()	{ return scoredSpec2; }
	protected BacktrackTable2D<T> getBacktrackTable()	{ return backtrackTable; }
	public AminoAcidSet getAASet()	{ return aaSet; }

	public float getScoreForDBScan(Annotation annotation)
	{
		int score1 = getScore1(annotation);
		int score2 = getScore2(annotation);
//		return score1+score2;
		float combinedScore = (float)(-Math.log(getSpectralProbability1(score1))-Math.log(getSpectralProbability2(score2)));
		return combinedScore;
	}
	
	public int getScore1(Annotation annotation) {
//		if(annotation == null)
//			return scoredSpec1.getScoreMinThreshold()-1;
//		if(srmScore1 != null)
//		{
//			Sequence<T> seq = graph.toCumulativeSequence(false, pep);
//			int score = 0;
//
//			for(int i=0; i<seq.size()-1; i++)
//			{
//				T srm = seq.get(i);
//				score += srmScore1.get(srm);
//			}
//			return score;
//		}
//		else
//			return scoredSpec1.getScore(pep, graph);
		return -1;
	}
	
	public int getScore2(Annotation annotation) {
//		if(pep == null)
//			return scoredSpec2.getScoreMinThreshold()-1;
//		if(srmScore2 != null)
//		{
//			Sequence<T> seq = graph.toCumulativeSequence(false, pep);
//			int score = 0;
//			for(int i=0; i<seq.size()-1; i++)
//			{
//				T srm = seq.get(i);
//				score += srmScore2.get(srm);
//			}
//			return score;
//		}
//		else
//			return scoredSpec2.getScore(pep, graph);
		return -1;
	}
	
	public int getScore1(Peptide pep)
	{
		return -1;
	}
	
	public int getScore2(Peptide pep)
	{
		return -1;
	}
	
	public int getMaxScore1()	{ return distribution.getMaxScore1(); }
	public int getMaxScore2()	{ return distribution.getMaxScore2(); }
	
	public int getMinScore1()	{ return distribution.getMinScore1(); }
	public int getMinScore2()	{ return distribution.getMinScore2(); }
	
	public float getSpectralProbability1(Peptide seq)
	{
		int score1 = getScore1(seq);
		return distribution.getSpectralProbability1(score1);
	}
	
	public float getSpectralProbability1(int score1)
	{
		return distribution.getSpectralProbability1(score1);
	}
	
	public float getSpectralProbability2(Peptide seq)
	{
		int score2 = getScore2(seq);
		return distribution.getSpectralProbability2(score2);
	}
	
	public float getSpectralProbability2(int score2)
	{
		return distribution.getSpectralProbability2(score2);
	}
	
	public float getSpectralProbabilitySumScores(Peptide seq) {
		int score1 = getScore1(seq);
		int score2 = getScore2(seq);
		return distribution.getSpectralProbabilitySumScores(score1, score2);
	}
	
	public float getNumEqualBetterPeptidesSumScores(Peptide seq)
	{
		int score1 = getScore1(seq);
		int score2 = getScore2(seq);
		return distribution.getNumEqualOrBetterPeptidesSumScores(score1, score2);
	}
	
	public float getSpectralProbability(Peptide seq) {
		int score1 = getScore1(seq);
		int score2 = getScore2(seq);
		return distribution.getSpectralProbability(score1, score2);
	}

	public float getProbabilityAt(int score1, int score2) {
		return distribution.getProbabilityAt(score1, score2);
	}
	
	public float getNumRecsAt(int score1, int score2) {
		return distribution.getNumRecsAt(score1, score2);
	}
	
	// score: inclusive
	public float getSpectralProbability(int score1, int score2) {
		return distribution.getSpectralProbability(score1, score2);
	}
	
	public float getProbBetterBoth(Peptide seq) {
		int score1 = getScore1(seq);
		int score2 = getScore2(seq);
		return distribution.getProbBetterBoth(score1, score2);
	}
	
	public float getProbBetterBoth(int score1, int score2) {
		return distribution.getProbBetterBoth(score1, score2);
	}
	
	public float getNumBetterBoth(Peptide seq) {
		int score1 = getScore1(seq);
		int score2 = getScore2(seq);
		return distribution.getNumBetterBoth(score1, score2);
	}
	
	public float getNumBetterBoth(int score1, int score2) {
		return distribution.getNumBetterBoth(score1, score2);
	}
	
	public float getNumEqualBetterPeptides(Peptide seq)
	{
		int score1 = getScore1(seq);
		int score2 = getScore2(seq);
		return distribution.getNumEqualOrBetterPeptides(score1, score2);
	}

	public float getNumEqualBetterPeptides1(Peptide seq)
	{
		return distribution.getNumEqualOrBetterPeptides1(getScore1(seq));
	}

	public float getNumEqualBetterPeptides2(Peptide seq)
	{
		return distribution.getNumEqualOrBetterPeptides2(getScore2(seq));
	}
	
	public float getNumEqualOrBetterPeptides(int score1, int score2)
	{
		return distribution.getNumEqualOrBetterPeptides(score1, score2);
	}

	// returns t where totalProb(t) > specProb && totalProb(t+1) <= specProb
	public int getThresholdScore(float specProb)
	{
		return -1;
	}
	
	/**
	 * Generate reconstructions with score "score" and have match with "sa" and put it in "reconstructions".
	 * @param score the score of reconstructions to be generated
	 * @param reconstructions a container where reconstructions will be stored
	 * @param sa suffix array that will filter reconstructions
	 * @return
	 */
	private void generateReconstructions(int score1, int score2, ArrayList<String> reconstructions, SuffixArray sa)
	{
		if(backtrackTable == null)
			return;
		for(T sink : graph.getSinkList())
			backtrackTable.getReconstructions(sink, score1, score2, "", reconstructions, sa);
	}

	public String getOneReconstruction(int score1, int score2)
	{
		if(backtrackTable == null)
			return null;
		return backtrackTable.getOneReconstruction(graph.getPMNode(), score1, score2, "");
	}
	
	public ArrayList<String> getReconstructions(int score1, int score2)
	{
		ArrayList<String> reconstructions = new ArrayList<String>();
		generateReconstructions(score1, score2, reconstructions, null);
		return reconstructions;
	}
	
	public ArrayList<String> getReconstructionsEqualOrAboveScore(int score1, int score2)
	{
		ArrayList<String> reconstructions = new ArrayList<String>();
		for(int t1=this.getMaxScore1()-1; t1>=score1; t1--)
			for(int t2=this.getMaxScore2()-1; t2>=score2; t2--)
				generateReconstructions(t1, t2, reconstructions, null);
		return reconstructions;
	}
	
	public ArrayList<String> getDictionary(float specProbThreshold)
	{
		return null;
	}

	public ArrayList<String> getReconstructions(float specProbThreshold, float numRecsThreshold, boolean isNumInclusive, SuffixArray sa)
	{
		return null;
	}
	
	public void computeGeneratingFunction()
	{
		// initialization of the source
		ScoreBound2D sourceBound = new ScoreBound2D(0,1,0,1);
		ScoreDist2D source = new ScoreDist2D(sourceBound);
		if(calcNumber)
			source.setNumber(0, 0, 1);
		if(calcProb)	
			source.setProb(0, 0, 1);
		
		fwdTable = new Hashtable<T, ScoreDist2D>();
		fwdTable.put(graph.getSource(), source);
		if(backtrack)
		{
			backtrackTable = new BacktrackTable2D<T>(graph, aaSet);
			BacktrackPointer2D sourcePointer = new BacktrackPointer2D(sourceBound, 0, 0);
			sourcePointer.setBacktrack(0, 0, 0);
			backtrackTable.put(graph.getSource(), sourcePointer);
		}
		
		// dynamic programming, source node (i=0) is excluded
		ArrayList<T> intermediateNodeList = graph.getIntermediateNodeList();
		// scores are ca
		
		srmScore1 = new Hashtable<T, Integer>();
		srmScore2 = new Hashtable<T, Integer>();
		
		for(int i=1; i<intermediateNodeList.size(); i++)
		{
			T curNode = intermediateNodeList.get(i);
			T srm = curNode;	// SRM
			T prm = graph.getComplementNode(srm);	// PRM
			int curScore1 = scoredSpec1.getNodeScore(prm, srm);
			srmScore1.put(curNode, curScore1);
			int curScore2 = scoredSpec2.getNodeScore(prm, srm);
			srmScore2.put(curNode, curScore2);
			
			setCurNode(curNode, curScore1, curScore2);
		}
		
		// process dest node
		int minScore1 = Integer.MAX_VALUE;
		int minScore2 = Integer.MAX_VALUE;
		int maxScore1 = Integer.MIN_VALUE;
		int maxScore2 = Integer.MIN_VALUE;
		
		for(T curNode : graph.getSinkList())
		{
			setCurNode(curNode, 0, 0);
			ScoreDist2D curDist = fwdTable.get(curNode);
			if(curDist == null)	// curNode is not connected from source (possible only when trypticOnly==true
				continue;
			if(curDist.getMinScore1() < minScore1)
				minScore1 = curDist.getMinScore1();
			if(curDist.getMinScore2() < minScore2)
				minScore2 = curDist.getMinScore2();
			if(curDist.getMaxScore1() > maxScore1)
				maxScore1 = curDist.getMaxScore1();
			if(curDist.getMaxScore2() > maxScore2)
				maxScore2 = curDist.getMaxScore2();
		}		

		assert(maxScore1 > minScore1);
		assert(maxScore2 > minScore2);
		
		ScoreDist2D finalDist = new ScoreDist2D(minScore1, maxScore1, minScore2, maxScore2); 
		for(T sinkNode : graph.getSinkList())
		{
			if(calcNumber)
				finalDist.addNumDist(fwdTable.get(sinkNode), 0, 0, 1);
			if(calcProb)
				finalDist.addProbDist(fwdTable.get(sinkNode), 0, 0, 1);
		}
		this.distribution = new ScoreDistMerged(finalDist);
	}
	
	private void setCurNode(T curNode, int curScore1, int curScore2)
	{
		int prevMinScore1 = Integer.MAX_VALUE;
		int prevMinScore2 = Integer.MAX_VALUE;
		int prevMaxScore1 = Integer.MIN_VALUE;
		int prevMaxScore2 = Integer.MIN_VALUE;
		// determine minScore and maxScore
		for(Edge<T> edge : graph.getEdges(curNode))
		{
			T prevNode = edge.getPrevNode();
			ScoreDist2D prevDist = fwdTable.get(prevNode);
			if(prevDist != null)
			{
				if(prevDist.getMaxScore1() > prevMaxScore1)
					prevMaxScore1 = prevDist.getMaxScore1();
				if(prevDist.getMaxScore2() > prevMaxScore2)
					prevMaxScore2 = prevDist.getMaxScore2();
				if(prevDist.getMinScore1() < prevMinScore1)
					prevMinScore1 = prevDist.getMinScore1();
				if(prevDist.getMinScore2() < prevMinScore2)
					prevMinScore2 = prevDist.getMinScore2();
			}
		}
		if(prevMinScore1 >= prevMaxScore1 || prevMinScore2 >= prevMaxScore2)
			return;
//		assert(prevMinScore1 < prevMaxScore1 || prevMinScore2 < prevMaxScore2);
		if(prevMaxScore1-prevMinScore1 > numScoreBinsPerNode)
			prevMinScore1 = prevMaxScore1-numScoreBinsPerNode;
		if(prevMaxScore2-prevMinScore2 > numScoreBinsPerNode)
			prevMinScore2 = prevMaxScore2-numScoreBinsPerNode;
		// compute score distribution
//		ScoreDist curDist = new ScoreDist(prevMinScore+curScore, prevMaxScore+curScore, calcNumber, calcProb); 
		ScoreBound2D curBound = new ScoreBound2D(prevMinScore1+curScore1, prevMaxScore1+curScore1, prevMinScore2+curScore2, prevMaxScore2+curScore2);
		ScoreDist2D curDist = new ScoreDist2D(curBound);
		BacktrackPointer2D backPointer = null;
		if(backtrack)
			backPointer = new BacktrackPointer2D(curBound, curScore1, curScore2);
		for(DeNovoGraph.Edge<T> edge : graph.getEdges(curNode))
		{
			T prevNode = edge.getPrevNode();
			int edgeIndex = edge.getEdgeIndex();
			ScoreDist2D prevDist = fwdTable.get(prevNode);
			if(prevDist != null)
			{
				if(calcNumber)
					curDist.addNumDist(prevDist, curScore1, curScore2, 1);
				if(calcProb)
					curDist.addProbDist(prevDist, curScore1, curScore2, edge.getEdgeProbability());
				if(backtrack)
				{
					BacktrackPointer2D prevPointer = backtrackTable.get(prevNode);
					backPointer.addBacktrackPointers(prevPointer, edgeIndex);
				}
			}
		}
		fwdTable.put(curNode, curDist);
		if(backtrack)
			backtrackTable.put(curNode, backPointer);	
	}
}
