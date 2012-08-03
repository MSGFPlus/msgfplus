package edu.ucsd.msjava.msgf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import edu.ucsd.msjava.msgf.DeNovoGraph.Edge;
import edu.ucsd.msjava.msutil.Matter;
import edu.ucsd.msjava.msutil.Sequence;

//TODO: implement it again
public class ProfileGF<T extends Matter> {

	private final GeneratingFunction<T> gf;
	
	public ProfileGF(GeneratingFunction<T> gf)
	{
		this.gf = gf;
	}
	
	private HashMap<T, ScoreDist> bwdTable = null;
	private double sizeDictionary = 0;
	private Profile<T> profile = null;
	
	public HashMap<T, ScoreDist> getBwdTable() {return bwdTable;}
	
	public Sequence<NominalMass> getGappedPeptideWithNominalMasses(float scoreAbove, float profileThreshold)
	{
		if(bwdTable == null)
			return null;
		ProfileGF<T> templateProfileGF = new ProfileGF<T>(this.gf);
		Profile<NominalMass> templateProf = templateProfileGF.computeProfileOfScoreAboveTop(scoreAbove).getSpectralProfile().toNominalMasses();
		
		Sequence<NominalMass> template = templateProf.getNodesWithProbEqualOrHigherThan(0.99999f);
		
		Sequence<NominalMass> mask = getSpectralProfile().toNominalMasses().getNodesWithProbEqualOrHigherThan(profileThreshold);
		
		Sequence<NominalMass> gappedPeptide = Sequence.getIntersection(template, mask);
		
		return gappedPeptide;
	}
	
	public Sequence<T> getGappedPeptide(float templateFraction, float specProb, float profileThreshold)
	{
		if(bwdTable == null)
			return null;
		ProfileGF<T> templateProfile = new ProfileGF<T>(this.gf);
		Sequence<T> template = templateProfile.computeProfileOfScoreAboveTop(templateFraction).getSpectralProfile().getNodesWithProbEqualOrHigherThan(0.99f);
		Sequence<T> mask = this.computeProfile(specProb).getSpectralProfile().getNodesWithProbEqualOrHigherThan(profileThreshold);
		Sequence<T> gappedPeptide = Sequence.getIntersection(template, mask);
		
		return gappedPeptide;
	}
	
	public Profile<T> getSpectralProfile()
	{
		if(profile != null)
			return profile;
		
		if(gf.getFwdTable() == null || bwdTable == null)
			return null;
		Profile<T> profile = new Profile<T>();
		
		for(T m :bwdTable.keySet())
		{
			ScoreDist fwdDist = gf.getFwdTable().get(m);
			ScoreDist bwdDist = bwdTable.get(m);
			if(fwdDist != null && bwdDist != null)
			{
				int minScore = bwdDist.getMinScore();
				int maxScore = bwdDist.getMaxScore();
				float sumNumbers = 0;
				for(int t = minScore; t<maxScore; t++){
					double mult = bwdDist.getNumberRecs(t);
					if(mult != 0)
						sumNumbers += fwdDist.getNumberRecs(t)*mult;
				}
				if(sumNumbers > 0)
					profile.add(new ProfilePeak<T>(m, sumNumbers/(float)sizeDictionary));
			}
		}
		
		Collections.sort(profile);
		
		this.profile = profile;
		
		return this.profile;
	}
	
	public ProfileGF<T> computeProfileOfScoreAboveTop(float fraction)
	{
		int thresholdScore = Math.round((gf.getMaxScore()-1)*fraction);
		return computeProfile(thresholdScore);
	}
	
	public ProfileGF<T> computeProfileOfTopScoringPeptides()
	{
		int thresholdScore = gf.getMaxScore()-1;
		return computeProfile(thresholdScore);
	}
	
	public ProfileGF<T> computeProfile(float specProb)
	{
		
		int thresholdScore = gf.getThresholdScore(specProb)+1;
		if(thresholdScore >= gf.getMaxScore())
			thresholdScore = gf.getMaxScore()-1;
//		else if(thresholdScore < gf.getMaxScore()-gf.getNumScoreBinsPerNode())
//			thresholdScore = gf.getMaxScore()-gf.getNumScoreBinsPerNode();
		
		return computeProfile(thresholdScore);
	}
	
	// thresholdScore: inclusive
	public ProfileGF<T> computeProfile(int thresholdScore)
	{
		sizeDictionary = gf.getNumEqualOrBetterPeptides(thresholdScore);
		
		// backward dynamic programming table
		HashMap<T, ScoreDist> bwdTable = new HashMap<T, ScoreDist>();
		
		ScoreDistFactory factory = new ScoreDistFactory(true, false);
		
		// initialization of the sink nodes
		ArrayList<T> sinkList = gf.getGraph().getSinkList();
		for(T curNode : sinkList)
		{
			ScoreDist sinkFwd = gf.getFwdTable().get(curNode);
			if(sinkFwd != null && sinkFwd.getMaxScore() > thresholdScore)
			{
                ScoreDist bwdDist = factory.getInstance(thresholdScore, sinkFwd.getMaxScore()); //**
				for(int t = thresholdScore; t<bwdDist.getMaxScore(); t++)
				{
					bwdDist.setNumber(t, 1);
				}
				bwdTable.put(curNode, bwdDist);
			}
		}
		
		// process intermediate nodes
		ArrayList<T> intermediateNodeList = gf.getGraph().getIntermediateNodeList();
		// setup score bounds of the backward table
		for(int i=intermediateNodeList.size()-1; i>0; i--)
		{
			T curNode = intermediateNodeList.get(i);
			ScoreDist fwdDist = gf.getFwdTable().get(curNode);
			if(fwdDist != null)
			{
				ScoreDist bwdDist = factory.getInstance(fwdDist.getMinScore(), fwdDist.getMaxScore());
				bwdTable.put(curNode, bwdDist);
			}
		}		
		
		// backward dynamic programming
		// sink nodes
		for(int i=sinkList.size()-1; i>=0; i--)
		{
			T curNode = sinkList.get(i);
			setBackwardNodes(curNode, bwdTable);			
		}
		// intermediate/source nodes
		for(int i=intermediateNodeList.size()-1; i>0; i--)
		{
			T curNode = intermediateNodeList.get(i);
			setBackwardNodes(curNode, bwdTable);			
		}
		
	
		this.bwdTable = bwdTable;
		
		return this;
	}	

	//TODO: replace getPreviousNode
	private void setBackwardNodes(T curNode, HashMap<T, ScoreDist> bwdTable)
	{
		ScoreDist curBwdDist = bwdTable.get(curNode);
		if(curBwdDist == null)
			return;
		BacktrackPointer pointer = gf.getBacktrackTable().get(curNode);
		int curNodeScore = pointer.getNodeScore();

		int bits = 0;
		ScoreDist[] prevBwdDists = new ScoreDist[gf.getGraph().getAASet().size()];
		
		for(int score=curBwdDist.getMaxScore()-1; score>=curBwdDist.getMinScore(); score--)
		{
			double numberRecs = curBwdDist.getNumberRecs(score);
			if(numberRecs == 0) continue;
			for(Edge<T> edge : gf.getGraph().getEdges(curNode))
			{
				int aaIndex = edge.getEdgeIndex();
				T prevNode = edge.getPrevNode();
				
				if((bits & (1 << aaIndex)) == 0)
				{
					bits |= (1 << aaIndex);
					prevBwdDists[aaIndex] = bwdTable.get(prevNode);
				}
				ScoreDist prevBwdDist = prevBwdDists[aaIndex];
				if(prevBwdDist != null)
					prevBwdDist.addNumber(score-curNodeScore, numberRecs);
			}
//			for(int aaIndex : pointer.getBacktrackAAIndexList(score))
//			{
//				if((bits & (1 << aaIndex)) == 0){
//					bits |= (1 << aaIndex);
//					T prevNode = gf.getGraph().getPreviousNode(curNode, gf.getGraph().getAASet().getAminoAcid(aaIndex));
//					prevBwdDists[aaIndex] =  bwdTable.get(prevNode);
//				}
//				T prevNode = gf.getGraph().getPreviousNode(curNode, gf.getAASet().getAminoAcid(aaIndex));
//				ScoreDist prevBwdDist = prevBwdDists[aaIndex];
//				if(prevBwdDist != null)
//					prevBwdDist.addNumber(score-curNodeScore, numberRecs);
//			}
		}
	}
}
