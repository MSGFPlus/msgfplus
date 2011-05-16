package msscorer;

import msgf.ScoredSpectrum;
import msgf.Tolerance;
import msutil.IonType;
import msutil.Matter;
import msutil.Peak;
import msutil.Spectrum;

//TODO: Implement CachedScoredSpectrum
public class NewScoredSpectrum<T extends Matter> implements ScoredSpectrum<T> {

	private Spectrum spec;
	private NewRankScorer scorer;
	private Tolerance mme;
	private IonType[][] ionTypes;	// segmentNum, ionType
	private final int charge;
	private final float parentMass;
	private IonType mainIon;
	private Partition partition;	// partition of the last segment
	private float probPeak;
	
//	private HashMap<T, Float> prefixScore = null;
//	private HashMap<T, Float> suffixScore = null;

	//	// used if scorer supports edgeScore
//	private HashMap<T, Float> nodeMass = null;
	
	public NewScoredSpectrum(Spectrum spec, NewRankScorer scorer)
	{
		this.spec = spec;
		this.scorer = scorer;
		
		this.charge = spec.getCharge();
		this.parentMass = spec.getParentMass();
		this.mme = scorer.mme;
		
		int numSegments = scorer.getNumSegments();
		ionTypes = new IonType[numSegments][];
		for(int seg=0; seg<numSegments; seg++)
			ionTypes[seg] = scorer.getIonTypes(charge, parentMass, seg);

		// filter precursor peaks
		for(PrecursorOffsetFrequency off : scorer.getPrecursorOFF(spec.getCharge()))
			spec.filterPrecursorPeaks(mme, off.getReducedCharge(), off.getOffset());
		spec.setRanksOfPeaks();
		
		// for edge scoring
		partition = scorer.getPartition(spec.getCharge(), spec.getParentMass(), scorer.getNumSegments()-1);
		mainIon = scorer.getMainIonType(partition);
		
		float approxNumBins = spec.getPeptideMass()/(scorer.getMME().getValue()*2);
		probPeak = spec.size()/approxNumBins;
	}
	
	@Override
	public int getNodeScore(T prm, T srm) 
	{
		float prefScore = getNodeScore(prm, true);
		float sufScore = getNodeScore(srm, false);
		return Math.round(prefScore+sufScore);	
	}	
	
	@Override
	public int getEdgeScore(T curNode, T prevNode, float theoMass) {
		if(!scorer.supportEdgeScores())
			return 0;
		
		int ionExistenceIndex = 0;
		float curNodeMass = getNodeMass(curNode);
		if(curNodeMass >= 0)
			ionExistenceIndex += 1;
		Float prevNodeMass = getNodeMass(prevNode);
		if(prevNodeMass >= 0)
			ionExistenceIndex += 2;
		
		float edgeScore = scorer.getIonExistenceScore(partition, ionExistenceIndex, probPeak);
		if(ionExistenceIndex == 3)
			edgeScore += scorer.getErrorScore(partition, curNodeMass-prevNodeMass-theoMass);
		return Math.round(edgeScore);
	}

	public NewRankScorer getScorer()
	{
		return scorer;
	}
	
	public Partition getPartition()
	{
		return partition;
	}
	
	public float getProbPeak()
	{
		return probPeak;
	}
	
	public IonType getMainIon()
	{
		return mainIon;
	}
	
	@Override
	public boolean getMainIonDirection() {
		return mainIon.isPrefixIon();
	}
	
	/**
	 * returns the corrected mass of the node based on the peak observed in the spectrum
	 * @param node
	 * @return	corrected mass of the node if peak exists, null -1
	 */
	public float getNodeMass(T node)
	{
		if(node.getNominalMass() == 0)
			return 0;
		float theoMass = mainIon.getMz(node.getMass());
		Peak p = spec.getPeakByMass(theoMass, scorer.getMME());
		if(p != null)
			return mainIon.getMass(p.getMz());
		else
			return -1;
	}
	
	public float getNodeScore(T node, boolean isPrefix)
	{
		float nodeMass = node.getMass();
		float score = 0;
		for(int segIndex=0; segIndex<scorer.getNumSegments(); segIndex++)
		{
			for(IonType ion : ionTypes[segIndex])
			{
				float theoMass;
				if(isPrefix)	// prefix
				{
					if(ion instanceof IonType.PrefixIon)
						theoMass = ion.getMz(nodeMass);
					else
						continue;
				}
				else
				{
					if(ion instanceof IonType.SuffixIon)
						theoMass = ion.getMz(nodeMass);
					else
						continue;
				}
				
				int segNum = scorer.getSegmentNum(theoMass, parentMass);
				if(segNum != segIndex)
					continue;
				
				Peak p = spec.getPeakByMass(theoMass, mme);
				Partition part = scorer.getPartition(charge, parentMass, segNum);
				
				if(p != null)	// peak exists
					score += scorer.getNodeScore(part, ion, p.getRank());
				else	// missing peak
					score += scorer.getMissingIonScore(part, ion);
			}			
		}
		
		return score;
	}

//	public void precomputeNodeScores(ArrayList<T> nodes)
//	{
//		if(prefixScore != null)
//			return;
//		prefixScore = new HashMap<T, Float>();
//		suffixScore = new HashMap<T, Float>();
//		
//		// intermediate nodes
//		for(T node : nodes)
//		{
//			float prefScore = getNodeScore(node, true);
//			float sufScore = getNodeScore(node, false);
//			prefixScore.put(node, prefScore);
//			suffixScore.put(node, sufScore);
//		}
//	}
//	
//	public void precomputeNodeMasses(ArrayList<T> nodes)
//	{
//		if(!scorer.supportEdgeScores())
//			return;
//		if(nodeMass != null)
//			return;
//		
//		nodeMass = new HashMap<T, Float>();
//		
//		int numNodesWithPeaks = 0;
//		for(T node : nodes)
//		{
//			float theoMass = mainIon.getMz(node.getMass());
//			Peak p = spec.getPeakByMass(theoMass, scorer.getMME());
//			if(p != null)
//			{
//				nodeMass.put(node, mainIon.getMass(p.getMz()));
//				numNodesWithPeaks++;
//			}
//		}
//		
//		this.probPeak = numNodesWithPeaks/(float)nodes.size();
//	}	
}
