package msgap;


import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import msgf.AminoAcidGraph;
import msgf.BacktrackPointer;
import msgf.DeNovoGraph;
import msgf.GeneratingFunction;
import msgf.NominalMass;
import msgf.Profile;
import msgf.ProfileGF;
import msgf.ProfilePeak;
import msgf.ScoreDist;
import msgf.ScoreDistFactory;
import msgf.ScoredSpectrum;
import msgf.DeNovoGraph.Edge;
import msutil.*;

public class GappedGeneratingFunction<T extends Matter> extends GeneratingFunction<T>{
	ArrayList<T> hubSet = new ArrayList<T>();
	//ScoreDist distribution = null;

	GapBacktrackTable<T> backtrackTable = null;	
	ProfileGF<T> profGF = null;
	// dynamic programming table
	ScoreDist[] scoreDistInHubs = null;
	// gap distribution table 
	GapFeatureTable<T> gapFeatureTable = null;
	ScoreDist distWellCleaved = null;
	
	public static HashMap<Integer, ArrayList<String>> getAACombTable(){
		return null;
	}
	
	public GappedGeneratingFunction(ScoredSpectrum<T> scoredSpec, DeNovoGraph<T> graph, Parameters par){
		super(scoredSpec, graph);
		this.enzyme(par.enzyme());
		
		if(par.allowNonEnzymaticCleavage()){
			GeneratingFunction<T> gf = new GeneratingFunction<T>(scoredSpec, graph);
			gf.enzyme(par.enzyme());
			gf.doNotBacktrack().doNotCalcNumber().computeGeneratingFunction();
			distWellCleaved = gf.getScoreDist();
			graph.allowNonEnzymaticCleavage();
		}
		
		this.computeGeneratingFunction();
		generateHubSet(par.specProbThresholdBeforeConsideringFlankingAAs(), par.numHubs());
		
	}
	
	public ScoreDist getDistWellCleaved(){ return distWellCleaved; }
	
	public ScoredSpectrum<T> getScoredSpectrum() { return super.getScoredSpectrum(); }
	/*
	private void initAAMassTable(){
		//
		if(aaMassTable == null){
			//System.out.println("init AA mass table..");
			aaMassTable = new HashMap<Integer, AminoAcid>();
			{
				for(AminoAcid aa : super.getAASet()){
					if(aa.getResidue() == 'L' || aa.getResidue() == 'Q') continue;
					aaMassTable.put(aa.getNominalMass(), aa);
				}
			}
		}
		if(aaCombTable == null){
			//System.out.println("init AA comb table..");
			aaCombTable = new HashMap<Integer, ArrayList<String>>();
			ArrayList<String> s = new ArrayList<String>();
			s.add("");
			aaCombTable.put(0, s);
			
			for(int i=50; i<=500; i++){
				ArrayList<String> v = new ArrayList<String>();
				for(AminoAcid aa : super.getAASet()){
					int prevmass = i - aa.getNominalMass();
					if(aaCombTable.containsKey(prevmass)){
						for(String t : aaCombTable.get(prevmass)){
							v.add(t+aa.getResidue());
						}
					}
				}
				
				if(!v.isEmpty()) aaCombTable.put(i, v);
			}
			//System.out.println("done init AA comb table..");
		}
	}*/
	
	/*public int[] doRecsContainPeptide(ArrayList<int[]> deNovo, Peptide pep){
		int[] correctGappedPeptide = fromPeptidetoGappedPeptide(pep);
		for(int[] rec : deNovo){
			if(GappedGeneratingFunction.isAcoveredbyB(correctGappedPeptide, rec))
				return rec;
		}
		return null;
	}*/

	
	public void setBackTrackPara(boolean orderByCoverage, int gappedPeptideLengthLimit, int dictionarySizeLimit){
		backtrackTable.setBackTrackPara(orderByCoverage, gappedPeptideLengthLimit, dictionarySizeLimit);
	}
	
	public void setBackTrackPara(boolean orderByCoverage, int gappedPeptideLengthLimit, int dictionarySizeLimit, int maxGapMass){
		backtrackTable.setBackTrackPara(orderByCoverage, gappedPeptideLengthLimit, dictionarySizeLimit, maxGapMass);
	}
	
	public void setRemoveRidundantGappedPeptides(boolean t){
		backtrackTable.setRemoveRedundantGappedPeptides(t);
	}
	
	
	
	/*
	public float getReducedDictionaryProb(int score, int delta){
		if(backtrackTable == null) return -1;
		if(score >= getMaxScore()) return -1;
		if(score < distribution.getMinScore()) score = distribution.getMinScore();
		
		backtrackTable.setCalReducedDictionaryProb(delta);
		for(T sink : getGraph().getSinkList()){
			backtrackTable.getReconstructionsEqualOrAboveScore(sink, score, null);
		}
		
		backtrackTable.offCalReducedDictionaryProb();
		return getSpectralProbability(score) - backtrackTable.getTotalShortGappedPeptideProb();
	}
	*/
	
	public float getGappedDictionaryProb(){return backtrackTable.getDictionaryProb();}
	
	/*
	public void generateHubSet(Profile<T>  profile, float intensityThreshold){
		for(T sink: graph.getSinkList())
			hubSet.add(sink);
		ArrayList<T> intermediateNodeList = graph.getIntermediateNodeList();
		for(int i=intermediateNodeList.size()-1; i>0; i--){
			T curNode = intermediateNodeList.get(i);
			int suffixMass = graph.getNominalPeptideMass() -  curNode.getNominalMass();
			for(ProfilePeak<T> p : profile){
				if((int)p.getNode().getNominalMass() == suffixMass){
					if(p.getProbability() >= intensityThreshold) hubSet.add(curNode);
					break;
				}
			}
		}
		hubSet.add(graph.getSource());
		generateIntHubSet();
	}
	
	public void generateHubSet(Profile<T> profile){
		HashSet<Short> possibleMasses = aaMassCombinationtable.get((short)558);
		//System.out.println(aaMassCombinationtable.get((short)558));
		for(T sink: graph.getSinkList())
			hubSet.add(sink);
		ArrayList<T> intermediateNodeList = graph.getIntermediateNodeList();
		float [] probabilities = new float[intermediateNodeList.size()+1];
		int index = 0;
		for(ProfilePeak<T> p : profile){
			probabilities[index++] = p.getProbability();
		}
		Arrays.sort(probabilities);
		
		//if(probabilities.length < noBlackNodes)
			//noBlackNodes = probabilities.length;
		
		
		for(int i=intermediateNodeList.size()-1; i>0; i--){
			T curNode = intermediateNodeList.get(i);
			if(profile.getNodesWithProbEqualOrHigherThan(probabilities[probabilities.length - 3]).contains(curNode)){
				hubSet.add(curNode);
				if(graph.getNominalPeptideMass() - curNode.getNominalMass() > 558){
					for(int j=i; j< intermediateNodeList.size(); j++){
						T tmpNode = intermediateNodeList.get(j);
						if(possibleMasses.contains((short)(tmpNode.getNominalMass() - curNode.getNominalMass())))
							if(!hubSet.contains(tmpNode)) hubSet.add(tmpNode);
					}
				}else{
					for(int j=i; j > 0; j--){
						T tmpNode = intermediateNodeList.get(j);
						if(possibleMasses.contains((short)(curNode.getNominalMass() -tmpNode.getNominalMass())))
							if(!hubSet.contains(tmpNode)) hubSet.add(tmpNode);
					}
				}
				break;
			}
		}
		hubSet.add(graph.getSource());
		generateIntHubSet();
	}
	*/
	
	
	private void generateHubSet(float specProb, int numHubs){
		
	//	System.out.println("1" + (System.nanoTime() - starttime));starttime = System.nanoTime();
		profGF = new ProfileGF<T>(this); // gf <- GeneratingFunction object
	
		int scoreThreshold = getThresholdScore(specProb);
		while(getNumEqualOrBetterPeptides(scoreThreshold) == Float.POSITIVE_INFINITY) scoreThreshold++;
	
		profGF.computeProfile(scoreThreshold); // generate profile
	//	System.out.println("2" + (System.nanoTime() - starttime));starttime = System.nanoTime();
		Profile<T> profile = profGF.getSpectralProfile(); // get the spectral profile
		
		numHubs -= getGraph().getSinkList().size();
		for(T sink: getGraph().getSinkList())
			hubSet.add(sink);
		
		ArrayList<Float> negatedProbabilities = new ArrayList<Float>();
		for(int i=0; i<numHubs;i++)
			negatedProbabilities.add(0f);
		
		for(ProfilePeak<T> p : profile){
			float prob = p.getProbability();
			if(-negatedProbabilities.get(numHubs-1) > prob) continue;
			
			int index = Collections.binarySearch(negatedProbabilities, -prob);
			
			if(index >= 0)
				negatedProbabilities.add(index, -prob);
			else
				negatedProbabilities.add(-index-1,-prob);
		}
		
		Sequence<T> tmp = profile.getNodesWithProbEqualOrHigherThan(-negatedProbabilities.get(numHubs-1));
		
		for(int i=tmp.size()-1;i>=0;i--){
			if(!hubSet.contains(tmp.get(i)))// && hubSet.size() <= numHubs)
				hubSet.add(tmp.get(i));
			
			if(hubSet.size() > numHubs) break;
		}
		
		if(!hubSet.contains(getGraph().getSource()))
			hubSet.add(getGraph().getSource());
	
	}
	

	/*
	public int[] fromPeptidetoGappedPeptide(Peptide p){
		int intSfMass = p.getNominalMass();
		
		int[] gp = new int[p.size()];
		ArrayList<Integer> massList = new ArrayList<Integer>();
		int gpsize = 0;
		
		for(T aa: hubSet)
			massList.add(aa.getNominalMass());
		
		for(int i=0; i<p.size(); i++){
			AminoAcid aa = p.get(i);
			if(massList.contains(intSfMass) && intSfMass > 0) gp[gpsize++] = intSfMass;
			intSfMass -= aa.getNominalMass(); 
		}
		
		int[] ret = new int[gpsize];
		for(int i=0; i<ret.length; i++){ ret[i] = gp[i];
		
		}

		return ret;
	}
	*/

	/*
	public void generateHubSet(float prmThreshold){
		for(T sink: graph.getSinkList())
			hubSet.add(sink);
		ArrayList<T> intermediateNodeList = graph.getIntermediateNodeList();
		for(int i=intermediateNodeList.size()-1; i>0; i--){
			T curNode = intermediateNodeList.get(i);
			int curScore = scoredSpec.getScore(curNode, graph.getPMNode());
			if(curScore >= prmThreshold)  hubSet.add(curNode); 
		}
		hubSet.add(graph.getSource());
		generateIntHubSet();
	}
	
	public void generateHubSet(int noBlackNodes){
		noBlackNodes-= 2; // sink + source
		for(T sink: graph.getSinkList())
			hubSet.add(sink);
		ArrayList<T> intermediateNodeList = graph.getIntermediateNodeList();
		int[] scores = new int[intermediateNodeList.size() - 1]; 
		for(int i=1; i<intermediateNodeList.size(); i++){
			T curNode = intermediateNodeList.get(i);
			scores[i-1] = scoredSpec.getScore(curNode, graph.getPMNode());
		}
		Arrays.sort(scores);
		int prmThreshold = Integer.MIN_VALUE, blackNodeCounter = 0;
		
		if(noBlackNodes <= intermediateNodeList.size() - 1) 
			prmThreshold = scores[intermediateNodeList.size() - 1 - noBlackNodes];
		
		for(int i=intermediateNodeList.size()-1; i>0; i--){
			T curNode = intermediateNodeList.get(i);
			int curScore = scoredSpec.getScore(curNode, graph.getPMNode());
			if((curScore >= prmThreshold) && (blackNodeCounter < noBlackNodes)){
				hubSet.add(curNode);
				blackNodeCounter++;
			}
		}
		hubSet.add(graph.getSource());
		generateIntHubSet();
	}
	*/
	
	
	public int getNumHubs(){return hubSet.size();}
	
	public ArrayList<T> getHubSet(){return hubSet;}
	

	/*public float[] getSuffixMassofGappedPeptide(short[] gappedPeptide){
		float[]  srm = new float[gappedPeptide.length]; 
		int j = 0, k=0;
		for(int i=0; i<hubSet.size(); i++){
			for(;j<gappedPeptide.length;j++){
				int mass = hubSet.get(i).getNominalMass();
				if(mass > gappedPeptide[j]) break;
				if(mass == gappedPeptide[j]) {
					srm[k++] = hubSet.get(i).getMass();
					break;
				}
			}
		}
		return srm;
	}*/
	
	public ArrayList<GappedReconstruction<T>> getGappedReconstructionsEqualOrAboveScore(int score)
	{
		return getGappedReconstructionsEqualOrAboveScore(score, null);
	}
	
	public ArrayList<GappedReconstruction<T>> getGappedReconstructionsEqualOrAboveScore(int score, GappedGeneratingFunction<T> other)
	{
		if(backtrackTable == null) return null;
		//if(score >= getMaxScore()) return null;
		if(score < getScoreDist().getMinScore()) score =  getScoreDist().getMinScore();
		
		ArrayList<GappedReconstruction<T>> reconstructions = new ArrayList<GappedReconstruction<T>>();
		
		if(other != null) backtrackTable.addRecsFrom(other.backtrackTable);
		
		for(T sink : getGraph().getSinkList())
			reconstructions.addAll(backtrackTable.getGappedReconstructionsEqualOrAboveScore(sink, score));
			
		return reconstructions;
	}
	
	
	/*public float getNumEqualOrBetterGappedPeptides(int score){
		if(!this.gapNumDistribution.isNumSet()) return -1;
		return gapNumDistribution.getNumEqualOrBetterPeptides(score);
	}*/
	
	
	
	/*public float getDictionarySize(float specProb)
	{
		return getNumEqualOrBetterGappedPeptides(getThresholdScore(specProb));
	}*/

	// modified by Sangtae, getNominalPeptideMass() is no longer available.
	public int getNominalParentMass(){
//		return getGraph().getNominalPeptideMass();
		return getGraph().getPMNode().getNominalMass();
	}
	
	public void computeGappedGeneratingFunction()
	{
	//	int scoreMinThreshold = getScoredSpectrum().getScoreMinThreshold();
		
		ScoreDistFactory factory = new ScoreDistFactory(false, calcProb());
		
		scoreDistInHubs = new ScoreDist[hubSet.size()];
		gapFeatureTable = new GapFeatureTable<T>(hubSet.size());

		// initialization of the source
        ScoreDist sourceDist = factory.getInstance(0, 1);

		if(calcProb())	
			sourceDist.setProb(0, 1);

		scoreDistInHubs[hubSet.indexOf(getGraph().getSource())] = sourceDist;
		
		//gap distributions should be 0 matrices. Do not need source distribution.
		if(backtrack())
		{
			backtrackTable = new GapBacktrackTable<T>(this, gapFeatureTable, scoreDistInHubs);
			GapBacktrackPointer<T> sourcePointer = new GapBacktrackPointer<T>(hubSet, 0, 1);
			sourcePointer.setBacktrack(0, 0,  getGraph().getSource());
			//System.out.println(hubSet.indexOf(getGraph().getSource()) + " "+ hubSet.size());
			backtrackTable.put(getGraph().getSource(), sourcePointer);
		}
		
		// dynamic programming, source node (i=0) is excluded
		//ArrayList<T> intermediateNodeList = getGraph().getIntermediateNodeList();
		for(ProfilePeak<T> p : profGF.getSpectralProfile())
		//for(int i=1; i<intermediateNodeList.size(); i++)
		{
			T curNode = p.getNode();
			if(getGraph().getSinkList().contains(curNode)) continue;
			//T curNode = intermediateNodeList.get(i);
			setCurNode(curNode, getSRMScore().get(curNode), factory);	
		}
		
		// process dest node
		//int minScore = Integer.MAX_VALUE;
		//int maxScore = scoreMinThreshold;
		for(T sinkNode : getGraph().getSinkList()) // no need to change as sink nodes are hub
		{
			setCurNode(sinkNode, 0, factory);
			/*ScoreDist curDist = scoreDistInHubs[hubSet.indexOf(sinkNode)];
			if(curDist == null)	// curNode is not connected from source (possible only when trypticOnly==true
				continue;
			if(curDist.getMinScore() < minScore)
				minScore = curDist.getMinScore();
			if(curDist.getMaxScore() > maxScore)
				maxScore = curDist.getMaxScore();*/
		}		

		//System.out.println(maxScore + "\t" + minScore);
		//System.out.println(super.getNumEqualOrBetterPeptides(super.getDistribution().getMaxScore()-1) + "\t"+ super.getSpectralProbability(super.getDistribution().getMaxScore()-1) + "\t" + super.getDistribution().getMaxScore() + "\t" + super.getDistribution().getMinScore());
		
		/*
		ScoreDist finalDist = factory.getInstance(minScore, maxScore);
		
		for(T sinkNode : getGraph().getSinkList())
		{
			int sinkNodeIndex = hubSet.indexOf(sinkNode);
			if(calcProb())
				finalDist.addProbDist(scoreDistInHubs[sinkNodeIndex], 0, 1);
		}
		*/
	//	this.distribution = finalDist;
	}
	
	private void setCurNode(T curNode, int curScore, ScoreDistFactory factory)
	{

		ArrayList<Edge<T>> hubConnectedEdges = new ArrayList<Edge<T>>();
		ArrayList<Edge<T>> nonHubConnectedEdges = new ArrayList<Edge<T>>();
		HashMap<Integer, int[]> connectingGapScoreRange = new HashMap<Integer, int[]>();
		
		// modified by Sangtae
//		for(T prevNode : getGraph().getPreviousNodes(curNode)){
		for(DeNovoGraph.Edge<T> edge : getGraph().getEdges(curNode))
		{
			T prevNode = edge.getPrevNode();
		
			if(hubSet.contains(prevNode)){     // if previous node = hub
				hubConnectedEdges.add(edge);		
			}else if(gapFeatureTable.containsNode(prevNode)){
				nonHubConnectedEdges.add(edge);	
			}	
		}

	
		if(nonHubConnectedEdges.isEmpty() && hubConnectedEdges.isEmpty()) return;

		//////
		// determine prev maxs mins of gap distributions and distribution
		//////
		
		for(int hubIndex=hubSet.size()-1; hubIndex>=0; hubIndex--){
			if(scoreDistInHubs[hubIndex] == null) continue;
			T hub = hubSet.get(hubIndex);
			if(curNode.getMass() <= hub.getMass()) break;
		
			int prevMaxGapScore = Integer.MIN_VALUE;
			int prevMinGapScore = Integer.MAX_VALUE;
			boolean hubConnected = false;
			
			for(Edge<T> e : hubConnectedEdges){
				if(e.getPrevNode().equals(hub)){
					hubConnected = true;
					break;
				}
			}
			
			if(hubConnected){
				// Number of edge between connected hub node ~ current = 1
				if(prevMaxGapScore < 1) prevMaxGapScore = 1;
				if(prevMinGapScore > 0) prevMinGapScore = 0;
			}
				
			for(Edge<T> e : nonHubConnectedEdges){ // from all the hubs ~ connected white
				T nonHub = e.getPrevNode();
				ScoreDist prevGapDist = gapFeatureTable.getDistBetween(hubIndex, nonHub);
				if(prevGapDist != null){
					int gapDistmax = prevGapDist.getMaxScore();
					int gapDistmin = prevGapDist.getMinScore();
					prevMaxGapScore = gapDistmax > prevMaxGapScore?  gapDistmax : prevMaxGapScore;
					prevMinGapScore = gapDistmin < prevMinGapScore? gapDistmin : prevMinGapScore;
				}
			}
			
			if(prevMinGapScore >= prevMaxGapScore) continue;
			
			int[] range = {prevMinGapScore, prevMaxGapScore};
			connectingGapScoreRange.put(hubIndex, range); // now we have gap range
		}
		
		if(connectingGapScoreRange.isEmpty()) return; // current node has no meaning

		//////
		//update Gap dist table
		//////
		
		for(int hubIndex : connectingGapScoreRange.keySet()){
			int[] range = connectingGapScoreRange.get(hubIndex);
			ScoreDist curGapDist = factory.getInstance(curScore + range[0], curScore + range[1]);
			ModifiedAAinGap modifiedAADist = null;
		
			if(ModifiedAAinGap.isAnyAAModified())
				modifiedAADist = new ModifiedAAinGap();
			
			T hub = hubSet.get(hubIndex);
				// Number of edge between connected hub node ~ current = 1
			
			Edge<T> he = null;
			
			for(Edge<T> e : hubConnectedEdges){
				if(e.getPrevNode().equals(hub)){
					he = e;
					break;
				}
			}
			
			if(he != null){
				if(calcProb()) curGapDist.addProb(curScore, he.getEdgeProbability());
				if(ModifiedAAinGap.isAnyAAModified()){
					for(AminoAcid aa : getAASet().getNominalMassAASet().getAAListWithNominalMass(curNode.getNominalMass()-hub.getNominalMass()))
						modifiedAADist.setModifiedAANumber(aa);
				}
			}	
			
			for(Edge<T> e : nonHubConnectedEdges){ // from all the hubs ~ connected white
				T nonHub = e.getPrevNode();
				
				ScoreDist prevGapDist = gapFeatureTable.getDistBetween(hubIndex, nonHub);
				if(prevGapDist == null) continue;
				
				if(calcProb())
					//Pi(t) += Pj(t-curScore) * prob(a.a)
					curGapDist.addProbDist(prevGapDist, curScore, e.getEdgeProbability());

				if(ModifiedAAinGap.isAnyAAModified()){
					ModifiedAAinGap prevModifiedAADist = gapFeatureTable.getModifiedAADistBetween(hubIndex, nonHub);
					for(AminoAcid aa : getAASet().getNominalMassAASet().getAAListWithNominalMass(curNode.getNominalMass()-nonHub.getNominalMass()))
						modifiedAADist.addModifiedAADist(prevModifiedAADist, aa);
				}
			}
			gapFeatureTable.putDist(hubIndex, curNode, curGapDist); 	
			if(ModifiedAAinGap.isAnyAAModified()) gapFeatureTable.putModifiedAADist(hubIndex, curNode, modifiedAADist);
		}
		
		//////
		//update distribution table if current node is hub
		//////
		
		if(hubSet.contains(curNode)){
			ScoreDist curDist = getFwdTable().get(curNode); //**
			GapBacktrackPointer<T> backPointer = null;
				
			if(backtrack())
				backPointer = new GapBacktrackPointer<T>(hubSet, curDist.getMinScore(), curDist.getMaxScore());
			
			for(int hubIndex : gapFeatureTable.getConnectedHubIndicesFrom(curNode)){ 	
				T hub = hubSet.get(hubIndex);
				
				ScoreDist gapDist = gapFeatureTable.getDistBetween(hubIndex, curNode);// curCol.get(hubNode);
			//	if(prevDist == null || gapDist == null)
			//		System.out.println("**");
			//	if(prevDist != null && gapDist != null){
					// from Gap distribution min ~ max score
					// scorediff should be increasing
					/*
					for(int scorediff = gapDist.getMinScore(); scorediff < gapDist.getMaxScore();scorediff++){
						num = gapDist.getNumberRecs(scorediff);
						if(num == 0f) continue;
						curDist.addNumDist(prevDist, scorediff, (int) num);
						if(calcProb()){
							prob = gapDist.getProbability(scorediff);
							//if(prob != 0.0)
							curDist.addProbDist(prevDist, scorediff, prob);	
						}
						//if(calcGapNum) curDist.addNumGapDist(prevDist, scorediff);
					}
					*/
				if(backtrack()) backPointer.addBacktrackPointers(backtrackTable.get(hub), hubIndex, gapDist.getMaxScore()-1);
					
			}
			
			scoreDistInHubs[hubSet.indexOf(curNode)] = curDist;// put current node and current distribution in the table
			if(backtrack())  backtrackTable.put(curNode, backPointer); // put backtrack pointer	
		}
		
	}	

	
	
	
	
	/*private void addBitSetsWithCardinality(int cardinality, int size, ArrayList<BitSet> bl){
		addBitSetsWithCardinality(cardinality, size, 0,  new BitSet(), bl);
	}
	
	private void addBitSetsWithCardinality(int cardinality, int size, int currentloc, BitSet b, ArrayList<BitSet> bl){
		if(b.cardinality() == cardinality){
			bl.add(b); return;
		}
		if(b.cardinality() > cardinality) return;
		if(size - currentloc  < cardinality - b.cardinality()) return;

		BitSet c = new BitSet();
		c.or(b);
		addBitSetsWithCardinality(cardinality, size, currentloc+1, c, bl);
		BitSet d = new BitSet();
		d.or(b);
		d.set(currentloc);
	
		addBitSetsWithCardinality(cardinality, size, currentloc+1, d , bl);
		
	}
	// input : AGPBP  output : AGPBP AGpBP AGpBp AGPBp
		
	private void addMatcheswithPTMs(String pep, HashSet<String> matches, int num_p_to_write){ // only p and P

		int numP = 0;
		int[] pIndex = new int[pep.length()];
		for(int i=pep.indexOf('.'); i<pep.lastIndexOf('.'); i++)
			if(pep.charAt(i) == 'P') pIndex[numP++] = i;
		
		//System.out.println(numP + "\t" + num_p_to_write);
		if(numP < num_p_to_write) num_p_to_write = numP;
		
		if(num_p_to_write == 0){
			matches.add(pep);
			return;
		}
		
		ArrayList<BitSet> bl = new ArrayList<BitSet>();
		
		addBitSetsWithCardinality(num_p_to_write, numP, bl);
		
		for(BitSet bs : bl){
			char[] p = pep.toCharArray();
			 for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i+1)) {
			     p[pIndex[i]] = 'p';
			 }
			 String toadd = new String(p);
			 matches.add(toadd);
		}

	}
	
	// opt later
	public HashSet<String> getMatchesWithPTMs(HashSet<String> matches){
		
		if(getAASet().getAminoAcid('p') != null){
			HashSet<String> matchesWithPTMs = new HashSet<String>();
			float mathDiff = getAASet().getAminoAcid('p').getMass() -  getAASet().getAminoAcid('P').getMass();
			for(String s : matches){
				ArrayList<AminoAcid> aaArray = new ArrayList<AminoAcid>(); // make a function later...
				
				for(int i=s.indexOf('.')+1; i<s.lastIndexOf('.');i++){
					char c = s.charAt(i);
					aaArray.add(getAASet().getAminoAcid(c));
				}
				
				Peptide pep = new Peptide(aaArray);
				int num_p_to_write =	Math.round((getGraph().getParentMass() - pep.getParentMass()) / mathDiff);
				
				if(num_p_to_write > 8) continue; // max 8 hydroxy prolines are allowd per a peptide
				
				addMatcheswithPTMs(s, matchesWithPTMs, num_p_to_write);
			}

			matches = matchesWithPTMs;
		}
		
		return matches;
	}

	public HashSet<String> getMatchedPeptidesAndProteinsviaGappedTags(ArrayList<int[]> deNovo, GappedHashtable table){
		
		HashSet<String> matches = new HashSet<String>();

		GappedTagSet gs = getGappedTagSet(deNovo);

		if(getAASet().getAminoAcid('p') != null){
			gs.addGappedTagswithPTMs();
		}
		
			
		for(GappedTag tag : gs){
			ArrayList<String> matchedPeptidesAndProteins = null;
			
			if(getAASet().getAminoAcid('p') != null){
			
				matchedPeptidesAndProteins = table.getMatchedPeptidesAndProteinswithPTMs(tag.getNTermMass(), table.getHashKey(tag), tag.getCTermMass());
				
			}else {
			
				matchedPeptidesAndProteins = table.getMatchedPeptidesAndProteins(tag.getNTermMass(), table.getHashKey(tag), tag.getCTermMass());
			}
			
			matches.addAll(matchedPeptidesAndProteins);
			
		}

	
		return matches;
	}
	
	public GappedTagSet getGappedTagSet(ArrayList<int[]> deNovo){
		return new GappedTagSet(deNovo, intHubSet, super.getAASet());
	}
	
	public float getSpecProbFixed(Peptide pep, String match){
		float specProb;
		
		if(getEnzyme().isCleaved(pep))
		{
			if(gfWellCleaved == null)
			{
				gfWellCleaved = new GeneratingFunction<T>(getScoredSpectrum(), getGraph()).enzyme(getEnzyme());
				gfWellCleaved.doNotBacktrack().doNotCalcNumber();
				
				gfWellCleaved.computeGeneratingFunction();
			}
			specProb = gfWellCleaved.getSpectralProbability(pep)*gfWellCleaved.getNeighboringAAProbability(match);
		}
		else
			specProb = super.getSpectralProbability(pep)*super.getNeighboringAAProbability(match);
		return specProb;
		
	}
	
	public boolean writeDBsearchResults(BufferedWriter out, float specProb, GappedHashtable table, Spectrum spec, boolean isSpecProbExclusive) throws IOException{
		return writeDBsearchResults(out, specProb, table, spec, isSpecProbExclusive, 0, true);
	}

	public boolean writeDBsearchResults(BufferedWriter out, float specProb, GappedHashtable table, Spectrum spec, boolean isSpecProbExclusive, int matchSocreThreshold, boolean useGappedTag) throws IOException{
		boolean isMatchFound = false; // found at least one match?
		backtrackTable.setRemoveRedundantGappedPeptides(false);
		
		ArrayList<GappedReconstruction<T>> gappedReconstructions = getGappedReconstructionsEqualOrAboveScore(Math.max(getThresholdScore(specProb), matchSocreThreshold));
		ArrayList<int[]> deNovo = new ArrayList<int[]>();
		
		for(GappedReconstruction<T> gappedReconstruction: gappedReconstructions)
			deNovo.add(gappedReconstruction.getGappedPeptide());
		
		HashSet<String> matches = new HashSet<String>();
		HashMap<String, ArrayList<String>> pep_pro = new HashMap<String, ArrayList<String>>();
		
		if(useGappedTag)
			matches = getMatchedPeptidesAndProteinsviaGappedTags(deNovo, table);
		else{
			matches = DBSearch.matchedPeptides(deNovo);
		}
		
		if(useGappedTag){
			for(String s: matches){
				String[] token = s.split("&&");
				ArrayList<String> pro = null;
				if(pep_pro.containsKey(token[0]))
					pro = pep_pro.get(token[0]);
				else 
					pro = new ArrayList<String>();
				
				pro.add(token[1]);
				pep_pro.put(token[0], pro);
			}
			matches.clear();
			matches.addAll(pep_pro.keySet());
		}
		
		matches = getMatchesWithPTMs(matches);
		
		int charge = spec.getPrecursorPeak().getCharge();
		
		for(String match : matches){
			ArrayList<AminoAcid> aaArray = new ArrayList<AminoAcid>();
			for(int i=match.indexOf('.')+1; i<match.lastIndexOf('.');i++){
				char c = match.charAt(i);
				aaArray.add(getAASet().getAminoAcid(c));
			}
			
			Peptide pep = new Peptide(aaArray);
			
			if(!getGraph().isNonEnzymaticCleavageAllowed() && !getEnzyme().isCleaved(pep)) continue;
			int psmScore = getScoredSpectrum().getScore(getGraph().toSequence(pep), getGraph().getPMNode());
			if(psmScore < Math.max(getThresholdScore(specProb), matchSocreThreshold))
				continue;	
			
			float specProbFixed =getSpecProbFixed(pep, match);
			
			if(isSpecProbExclusive && specProbFixed > specProb) continue;
	
			System.out.println(spec.getScanNum() + "\t" + spec.getAnnotationStr() + " and " + pep + "\t" + psmScore);
			if(useGappedTag){
				for(String pro: pep_pro.get(match))
					out.write(spec.getTitle()+"\t"+spec.getScanNum()+"\t"+spec.getPrecursorPeak().getMz()+"\t"+charge+"\t" + match + "\t" + pro + "\t" + psmScore + "\t" + specProbFixed + "\n");
			}else
				out.write(spec.getTitle()+"\t"+spec.getScanNum()+"\t"+spec.getPrecursorPeak().getMz()+"\t"+charge+"\t" + match + "\t" + psmScore + "\t" + specProbFixed + "\n");
			isMatchFound = true;
		}
		
		return isMatchFound;
	}*/

	
	
    static ArrayList<GappedReconstruction<NominalMass>> getGappedDictionary(ArrayList<GappedGeneratingFunction<NominalMass>> gappedGFList, int thresholdScore){
    	ArrayList<GappedReconstruction<NominalMass>> gappedReconstructions = null;
		GappedGeneratingFunction<NominalMass> other = null;
		for(GappedGeneratingFunction<NominalMass> gap : gappedGFList){
			gappedReconstructions = gap.getGappedReconstructionsEqualOrAboveScore(thresholdScore, other);	
			//if(gappedReconstructions.size() == 0) System.out.println("++" + gappedReconstructions.size());
			other = gap;
		}
		
		for(int index = 0; index < gappedReconstructions.size(); index++){
			GappedReconstruction<NominalMass> gappedReconstruction = gappedReconstructions.get(index);
			for(GappedReconstruction<NominalMass> modifiedGappedReconstruction : gappedReconstruction.getModifiedGappedReconstructions()){
				gappedReconstructions.add(index+1, modifiedGappedReconstruction);
				index++;
			}
		}
		
		//System.out.println(other.getDictionaryProb()+"\t"+other.getSpectralProbability(other.getThresholdScore(specProb)));
		
		return gappedReconstructions;
    }
    
    ArrayList<GappedReconstruction<T>> generateGappedDictionary(int thresholdScore){
    	ArrayList<GappedReconstruction<T>> gappedReconstructions = null;
    	
    	gappedReconstructions = getGappedReconstructionsEqualOrAboveScore(thresholdScore);
	
		for(int index = 0; index < gappedReconstructions.size(); index++){
			GappedReconstruction<T> gappedReconstruction = gappedReconstructions.get(index);
			for(GappedReconstruction<T> modifiedGappedReconstruction : gappedReconstruction.getModifiedGappedReconstructions()){
				gappedReconstructions.add(index+1, modifiedGappedReconstruction);
				index++;
			}
		}
			
		return gappedReconstructions;
    }
    
    
    
}
