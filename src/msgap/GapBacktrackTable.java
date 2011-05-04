package msgap;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

import msgf.ScoreDist;
import msutil.Matter;

public class GapBacktrackTable<T extends Matter> extends HashMap<T, GapBacktrackPointer<T>> {
	
	private class GappedReconstructionComparator implements Comparator<Object>{
		@SuppressWarnings("unchecked")
		public int compare(Object gp1, Object gp2){
			float prob1 = ((GappedReconstruction<T>) gp1).getProb();
			float prob2 = ((GappedReconstruction<T>) gp2).getProb();
			if(prob1 > prob2)
				return -1;
			else if(prob1 < prob2)
				return 1;
			else
				return 0;
		}
	}

	private static final long serialVersionUID = 1L;
	private boolean orderByCoverage = true;
	private int delta = 5;//0;
	private int dictionarySizeLimit = 25;//Integer.MAX_VALUE; 
	ArrayList<GappedReconstruction<T>> gappedReconstructions = null;
	private float dictionaryProbability = 0;
	private GapFeatureTable<T> gapFeatureTable = null;
	private ArrayList<T> hubSet = null;
	private GappedGeneratingFunction<T> gap = null;
	private boolean removeRedundantGappedPeptides = false;
	private ScoreDist[] scoreDistInHubs = null;
	private  ScoreDist[][] prevGapDists = null;
	private int maxGapMass = Integer.MAX_VALUE;
	
	public GapBacktrackTable(GappedGeneratingFunction<T> gap, GapFeatureTable<T> gapDistTable, ScoreDist[] scoreDistInHubs)
	{
		this.gapFeatureTable = gapDistTable;
		this.scoreDistInHubs = scoreDistInHubs;
		this.gap = gap;
		this.hubSet = gap.getHubSet();
	}
	
	public void setRemoveRedundantGappedPeptides(boolean t) {removeRedundantGappedPeptides = t;}
	
	
	public void setBackTrackPara(boolean orderByCoverage, int gappedPeptideLengthLimit, int dictionarySizeLimit){
		setMinGappedPeptideLength(gappedPeptideLengthLimit);
		setOrderByCoverage(orderByCoverage);	
		setDictionarySizeLimit(dictionarySizeLimit);
	}
	
	public void setBackTrackPara(boolean orderByCoverage, int gappedPeptideLengthLimit, int dictionarySizeLimit, int maxGapMass){
		this.maxGapMass = maxGapMass;
		setBackTrackPara(orderByCoverage, gappedPeptideLengthLimit, dictionarySizeLimit);
	}
	
	public void setDictionarySizeLimit(int dictionarySizeLimit){ this.dictionarySizeLimit = dictionarySizeLimit;	}
	
	public void setOrderByCoverage(boolean tf){ this.orderByCoverage = tf; }
	
	public void setMinGappedPeptideLength(int delta){ this.delta = delta; }
		
	public float getDictionaryProb(){ return dictionaryProbability; }
		
	
	public void addProbToGappedReconstruction(int index, float prob){
		GappedReconstruction<T> gappedReconstruction = gappedReconstructions.get(index);
		float addedProb = gappedReconstruction.getProb() + prob;
		removeGappedReconstruction(index);
		gappedReconstruction.setProb(addedProb);
		addGappedReconstruction(gappedReconstruction, getGappedReconstructionRank(gappedReconstruction));
	}
	
	
	public void addGappedReconstruction(GappedReconstruction<T> gappedReconstruction){
		addGappedReconstruction(gappedReconstruction, gappedReconstructions.size());
	}
	
	public void addGappedReconstruction(GappedReconstruction<T> gappedReconstruction, int rank){
		dictionaryProbability += gappedReconstruction.getProb();
		gappedReconstructions.add(rank, gappedReconstruction);
	}
	
	public void removeLastGappedReconstruction(){ removeGappedReconstruction(gappedReconstructions.size() - 1); }
	
	public void removeGappedReconstruction(int index){
		dictionaryProbability -= gappedReconstructions.get(index).getProb();
		gappedReconstructions.remove(index);
	}
	
	
	private float calGappedPeptideProb(float[] prevSuffixProbs, float[] suffixProbs, int prevMinSuffixScore, int minSuffixScore, ScoreDist prevDist , ScoreDist gapDist ){
		float prevTotalProb = 0;
		//int maxSuffixScore = minSuffixScore + suffixProbs.length;
		for(int prevscore = prevMinSuffixScore; prevscore < prevDist.getMaxScore(); prevscore++){
			float prevScoreForwardProb = prevDist.getProbability(prevscore);
			if(prevScoreForwardProb == 0f) continue;
			float prevScoreProb = 0.0f;
			
			for(int gapScore =  Math.max(minSuffixScore-prevscore, gapDist.getMinScore()); gapScore < Math.min(gapDist.getMaxScore(), suffixProbs.length + minSuffixScore-prevscore); gapScore++){
				float gapProb =  gapDist.getProbability(gapScore);
				if(gapProb == 0f) continue;
				int suffixProbsIndex = prevscore + gapScore - minSuffixScore;
				
				if(suffixProbs[suffixProbsIndex] > 0)
					prevScoreProb += gapProb * suffixProbs[suffixProbsIndex];

					//prevScoreProb += gapDist.getNumberRecs(scorediff) * scoreProbs[score - minScore]; // to count the number of peptides that induce a gapped peptide
			}
			prevSuffixProbs[prevscore - prevMinSuffixScore] = prevScoreProb;
			// multiply by prevDist.getProbability(prevscore) to upper bound real prevTotalProb
			prevTotalProb += prevScoreProb * prevScoreForwardProb;
		}
		return prevTotalProb;
	}
	
	private int getGappedReconstructionRank(GappedReconstruction<T> gappedReconstruction){

		int index = Collections.binarySearch(gappedReconstructions, gappedReconstruction, new GappedReconstructionComparator());
		if(index < 0) return -index-1;
		else return index;
		
	}
	
	private int[] getNewAppendedHubIndices(int[] hubIndices, int hubIndex){
		int[] appendedHubIndices = new int[hubIndices.length + 1];
		for(int k=0; k<hubIndices.length; k++) appendedHubIndices[k] = hubIndices[k];
		appendedHubIndices[appendedHubIndices.length-1] = hubIndex;
		return appendedHubIndices;
	}
	
	public void backtrack(int hubIndex, int score, int[] hubIndices, float[] suffixProbs, int minSuffixScore, float totalSuffixProb)
	{
		
		GapBacktrackPointer<T> pointer = this.get(hubSet.get(hubIndex));
		if(pointer == null) return;
		if(score >= pointer.getMaxScore()) return;
		//if(calTotalShortGappedPeptideProb){ if(prefix.length >= delta) return; }
		if(!orderByCoverage) 
			if(gappedReconstructions.size() >= dictionarySizeLimit) return;
			
		
		if(hubIndex == hubSet.size() - 1){ // source
			if(hubIndices.length >= delta)
			{  
				// test		
				hubIndices = getNewAppendedHubIndices(hubIndices, hubIndex);
				GappedReconstruction<T> gappedReconstruction = new GappedReconstruction<T>(gap, hubIndices);
				
				if(maxGapMass < Integer.MAX_VALUE && gappedReconstruction.getMaxGapmass() > maxGapMass) return;
				
				//getMaxGapmass
				if(orderByCoverage){
					if(removeRedundantGappedPeptides){
						for(int i=0; i<gappedReconstructions.size(); i++){
							GappedReconstruction<T> other = gappedReconstructions.get(i);
							if(gappedReconstruction.isRedundantWRTOtherGappedReconstruction(other)){
								addProbToGappedReconstruction(i, totalSuffixProb); // prefix redundant 
								return;
							}						
							else if(other.isRedundantWRTOtherGappedReconstruction(gappedReconstruction)){
								totalSuffixProb += other.getProb();
								removeGappedReconstruction(i--);
							}
						}				
					}
					
					gappedReconstruction.setProb(totalSuffixProb);
					
					int rank = getGappedReconstructionRank(gappedReconstruction);
					if(rank >= dictionarySizeLimit) return;
					
					addGappedReconstruction(gappedReconstruction, rank);
					while(gappedReconstructions.size() > dictionarySizeLimit)
						removeLastGappedReconstruction();
				}else{
					gappedReconstruction.setProb(totalSuffixProb);
					addGappedReconstruction(gappedReconstruction);
				}
			}
			//}else if(prefix.length < delta) totalShortGappedPeptideProb +=totalSuffixProb;
			return;
		}
			
		
		float minGappedReconstructionProb = 0;
		int[] nextHubIndices = null;
		
		if(gappedReconstructions.size() >=  dictionarySizeLimit)
			minGappedReconstructionProb = gappedReconstructions.get(dictionarySizeLimit-1).getProb();
		
		for(int connectedHubIndex:pointer.getConnectedHubIndices(score))
		{	
			ScoreDist prevDist = scoreDistInHubs[connectedHubIndex];//dpTable.get(hub); 
			ScoreDist gapDist = prevGapDists[connectedHubIndex][hubIndex];//gapDistTable.getDistBetween(hub, curNode);
			int prevMinSuffixScore = Math.max(prevDist.getMinScore(), minSuffixScore - gapDist.getMaxScore() + 1);
			float[] prevSuffixProbs = new float [prevDist.getMaxScore() - prevMinSuffixScore];
			
			float prevTotalSuffixProb = calGappedPeptideProb(prevSuffixProbs, suffixProbs, prevMinSuffixScore, minSuffixScore, prevDist, gapDist);
			if(orderByCoverage) 
				if(prevTotalSuffixProb <= minGappedReconstructionProb)
					continue;
				
			if(nextHubIndices == null)
				nextHubIndices = getNewAppendedHubIndices(hubIndices, hubIndex);
			
			backtrack(connectedHubIndex, score - pointer.getScoreDiff(connectedHubIndex), nextHubIndices, prevSuffixProbs, prevMinSuffixScore, prevTotalSuffixProb);
		}
		return;		
	}
	
	public void addRecsFrom(GapBacktrackTable<T> other){
		if(other == null || other.gappedReconstructions == null) return;
		gappedReconstructions = new ArrayList<GappedReconstruction<T>>();
		gappedReconstructions.addAll(other.gappedReconstructions);
		dictionaryProbability += other.getDictionaryProb();
	}
	
	
	public ArrayList<GappedReconstruction<T>> getGappedReconstructionsEqualOrAboveScore(T sinkNode, int thresholdScore){
		
		if(gappedReconstructions == null) gappedReconstructions = new ArrayList<GappedReconstruction<T>>();
		
		if(prevGapDists == null){
			prevGapDists = new ScoreDist[hubSet.size()][hubSet.size()];
			
			for(int i=0;i<hubSet.size(); i++){
				for(int k=i+1;k<hubSet.size(); k++){
					prevGapDists[k][i] = gapFeatureTable.getDistBetween(k, hubSet.get(i));
				}
			}
		}
		// pre processing to calculate the probabilities of GPs
		int[] scoreBound = new int[2];
		ScoreDist sinkDist = scoreDistInHubs[hubSet.indexOf(sinkNode)];
		
		scoreBound[0] = thresholdScore; scoreBound[1] =  sinkDist.getMaxScore();
		if(scoreBound[1] < scoreBound[0]) return gappedReconstructions;
		float [] suffixProbs = new float [scoreBound[1] - scoreBound[0]]; 
		
	    for(int score = scoreBound[1] -1; score >= scoreBound[0]; score--)
			if(sinkDist.getNumberRecs(score) != 0.0)	suffixProbs[score - scoreBound[0]] = 1.0f;	
		
		// backtracking
	    for(int score = scoreBound[1] -1; score >= scoreBound[0]; score--){
			if(sinkDist.getNumberRecs(score) != 0.0)
				backtrack(hubSet.indexOf(sinkNode),  score, new int[0], suffixProbs, scoreBound[0], 1.0f);
		}
	    
	
	    //System.out.println(gappedReconstructions.size());
	    
	   /* for(int i=0; i<hubSet.size()-1;i++){
	    	if(gap.gapFeatureTable.getModifiedAADistBetween(i+1, hubSet.get(i)) != null){
	    	System.out.print(hubSet.get(i).getNominalMass()-hubSet.get(i+1).getNominalMass()+"\t" + gap.getAASet().getAminoAcid('c')+ "\t" + gap.getAASet().getAminoAcid('C')+"\t" + gap.getAASet().getAminoAcid('p')+ "\t" + gap.getAASet().getAminoAcid('P'));
	    	System.out.print(gap.gapFeatureTable.getModifiedAADistBetween(i+1, hubSet.get(i)).getModifiedAANumberMap().get(gap.getAASet().getAminoAcid('p')));
	    	System.out.print(gap.gapFeatureTable.getModifiedAADistBetween(i+1, hubSet.get(i)).getModifiedAANumberMap().get(gap.getAASet().getAminoAcid('c')));
	    	System.out.println("\n\t* " +gap.gapFeatureTable.getModifiedAADistBetween(i+1, hubSet.get(i)).getAllPossibleMassDifferences());
	    
	    	}
	    }*/
	    
	    ArrayList<Integer> toErase = new ArrayList<Integer>();
	    for(int i=0; i<gappedReconstructions.size(); i++){
	    	GappedReconstruction<T> gappedReconstruction = gappedReconstructions.get(i);
		    for(int j=0; j<gappedReconstructions.size(); j++){
		    	if(i==j) continue;
				GappedReconstruction<T> other = gappedReconstructions.get(j);
				if(gappedReconstruction.isRedundantWRTOtherGappedReconstruction(other)){
					toErase.add(i);
					break;
				}						
			}	//removeGappedReconstruction
	    }
	
	    Collections.sort(toErase);
	    for(int i=toErase.size()-1;i>=0;i--){
	    	removeGappedReconstruction(toErase.get(i));
	    }
	    
	    
	    return gappedReconstructions;
	}
}
