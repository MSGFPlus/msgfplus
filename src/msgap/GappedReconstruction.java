package msgap;

import java.util.ArrayList;
import java.util.HashMap;

import msutil.AminoAcid;
import msutil.Matter;
import msutil.Peptide;

public class GappedReconstruction <T extends Matter>{
	private GappedGeneratingFunction<T> gap = null;
	private int[] suffixMasses = null;
	private int[] hubIndices = null;
	private float prob;	
	private ArrayList<T> hubSet = null;
	private GappedReconstruction <T> gappedReconstructionBeforeModification = null;
	private boolean containsModification = false;

	GappedReconstruction(GappedGeneratingFunction<T> gap, int[] hubIndices){
		this.gap = gap;
		hubSet = gap.getHubSet();
		this.hubIndices = hubIndices;
		updatePrefixMassesFromHubIndices();
	}
	
	private GappedReconstruction(GappedGeneratingFunction<T> gap, GappedReconstruction<T> gappedReconstructionBeforeModification, ArrayList<Integer> gapMasses){
		this.gap = gap;
		hubSet = gap.getHubSet();
		hubIndices = gappedReconstructionBeforeModification.hubIndices;
		
		containsModification = true;
		this.gappedReconstructionBeforeModification = gappedReconstructionBeforeModification;
		
		suffixMasses = new int[gapMasses.size()];
		for(int i = gapMasses.size()-1;i>=0;i--){
			suffixMasses[i] = gapMasses.get(i) + (i<gapMasses.size()-1? suffixMasses[i+1] : 0);
		}
	}
	
	
	private void subGetModifiedGappedReconstructions(ArrayList<Integer> gapMasses, int modNum, ArrayList<ArrayList<Integer>> modifiedGappedReconstructionsRepresentedByGapMasses){
		int gapIndex = gapMasses.size();
		
		if(gapIndex >= suffixMasses.length){
			if(modNum>0){
				modifiedGappedReconstructionsRepresentedByGapMasses.add(gapMasses);
			}
			return;
		}
		
		ModifiedAAinGap md = gap.gapFeatureTable.getModifiedAADistBetween(hubIndices[gapIndex+1], hubSet.get(hubIndices[gapIndex]));
			//c,JJ,C5H8N1O3  p,HydroxyProline,C5H7N1O2

		if(md == null || modNum >= ModifiedAAinGap.getMaxModNum()){
			int nextGapMass = suffixMasses[gapIndex] - (gapIndex+1 < suffixMasses.length ? suffixMasses[gapIndex+1] : 0);
			ArrayList<Integer> nextGapMasses = new ArrayList<Integer>(gapMasses);
			nextGapMasses.add(nextGapMass);
			subGetModifiedGappedReconstructions(nextGapMasses, modNum, modifiedGappedReconstructionsRepresentedByGapMasses);
		}else{
			//if(!md.getAllPossibleMassDifferences().equals( md.getAllPossibleMassDifferences2())){
			//	System.out.print(prefixMasses[gapIndex] - (gapIndex+1 < prefixMasses.length ? prefixMasses[gapIndex+1] : 0) + " " + md.getAllPossibleMassDifferences());
			//	System.out.println(" new" + md.getAllPossibleMassDifferences2());
			//}
				
			for(int massDiff : md.getAllPossibleMassDifferences().keySet()){
				int nextGapMass = massDiff + suffixMasses[gapIndex] - (gapIndex+1 < suffixMasses.length ? suffixMasses[gapIndex+1] : 0);
				ArrayList<Integer> nextGapMasses = new ArrayList<Integer>(gapMasses);
				nextGapMasses.add(nextGapMass);
				subGetModifiedGappedReconstructions(nextGapMasses, modNum + md.getAllPossibleMassDifferences().get(massDiff), modifiedGappedReconstructionsRepresentedByGapMasses);
			}
		}
	}
	
	public ArrayList<GappedReconstruction<T>> getModifiedGappedReconstructions(){
		ArrayList<GappedReconstruction<T>> ret = new ArrayList<GappedReconstruction<T>>();
		ArrayList<ArrayList<Integer>> modifiedGappedReconstructionsRepresentedByGapMasses = new ArrayList<ArrayList<Integer>>();
		
		subGetModifiedGappedReconstructions(new ArrayList<Integer>(), 0, modifiedGappedReconstructionsRepresentedByGapMasses);
		
		for(ArrayList<Integer> gapMasses : modifiedGappedReconstructionsRepresentedByGapMasses){
			ret.add(new GappedReconstruction<T>(gap, this, gapMasses));
		}
	
		return ret;
	}
	
	float getProb(){ 
		if(containsModification)
			return gappedReconstructionBeforeModification.getProb();
		else return prob; 	
	}
	
	void addProb(float prob){
		this.prob += prob;
	}
	
	void setProb(float prob){
		this.prob = prob;
	}
	
	public int[] getGappedPeptide(){return suffixMasses;}
	
	void updatePrefixMassesFromHubIndices(){
		suffixMasses = new int[hubIndices.length-1];
		for(int i=0; i<hubIndices.length-1; i++)
			suffixMasses[i] = hubSet.get(hubIndices[i]).getNominalMass();
	}
	
	boolean isRedundantWRTOtherGappedReconstruction(GappedReconstruction<T> other){
		int[] b = other.suffixMasses;
		
		boolean ret = true;
		int a_index = 0;
		for(int i=0; i<b.length; i++){
			for(int j=a_index; j<suffixMasses.length; j++){
				if(b[i] == suffixMasses[j]){
					ret = true;
					a_index = j;
					break;
				}
				ret = false;
			}
			if(ret == false) break;			
		}
		return ret;
	}
	
	
	
	public ArrayList<Integer> getGapMassRepresentation(){
		ArrayList<Integer> gp = new ArrayList<Integer>();
		for(int i=0;i<suffixMasses.length-1;i++)
			gp.add(suffixMasses[i] - suffixMasses[i+1]);	
		gp.add((int)suffixMasses[suffixMasses.length-1]);
	
		return gp;
	}
	
	public int getLength() {
		return getGapMassRepresentation().size();
	}
	
	public ArrayList<Integer> getSuffixMassRepresentation(){
		ArrayList<Integer> gp = new ArrayList<Integer>();
		for(int m : suffixMasses){
			if(m>0)gp.add(m);
		}
		return gp;
	}
	
	public ArrayList<Integer> getPrefixMassRepresentation(){
		ArrayList<Integer> gp = new ArrayList<Integer>();
		int pm = suffixMasses[0];
		for(int m : suffixMasses){
			if(pm-m>0)gp.add(pm - m);
		}
		return gp;
	}
	
	public float getMaxGapmass(){
		ArrayList<Integer> gapMasses = getGapMassRepresentation();
		int max = 0;
		
		for(int m : gapMasses)
			max = Math.max(max, m);
		
		return max;
	}
	
	
	public String toString(){
		String ret = "REC:";
		
		for(int m : suffixMasses){
			ret+=" " + m;
		}
		
		return ret+"\n";
	}
	
	
	
}
