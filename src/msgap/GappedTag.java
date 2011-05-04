package msgap;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Hashtable;

import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Peptide;

public class GappedTag {
	private int length;
	private int numGap;
	private int maxGapSize = 0;
	private int maxGapSizeLimit = 3000;
	private int nTermMass, cTermMass;
	private float filteringEff = 1;
	char[] aaResidues = null;
	int[] gaps = null;
	private BitSet isAminoAcid = null;
	
	static AminoAcidSet aaSet = null;//AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
	static Hashtable<Integer, ArrayList<AminoAcid>> aaMassTable = null;
	static ArrayList<Integer> intHubSet = new ArrayList<Integer>();
	static void setIntHubSet(ArrayList<Integer> inthubset) {intHubSet = inthubset;}
	static Hashtable<Integer, Float> gapFilteringEffTable = new Hashtable<Integer, Float>();
	
	private void initAAMassTable(AminoAcidSet a){
		if(aaSet == null){ 
			aaSet = a;
			aaMassTable = new Hashtable<Integer, ArrayList<AminoAcid>>();
			for(AminoAcid aa : aaSet){
				if(aa.getResidueStr() == 'I') continue;
				
				short key = (short)aa.getNominalMass();
				ArrayList<AminoAcid> value = null;
				
				if(aaMassTable.containsKey(key))
					value = aaMassTable.get(key);
				else 
					value = new ArrayList<AminoAcid>();
				value.add(aa);
				
				aaMassTable.put(aa.getNominalMass(), value);
			}//aaSet.getProbability(aa)
			gapFilteringEffTable.put(0, 1.0f);
			
			for(int m=50; m<=maxGapSizeLimit;m++){
				float f = 0;
				for(int d: aaMassTable.keySet()){
					if(gapFilteringEffTable.containsKey(m-d)){
						for(AminoAcid aa : aaMassTable.get(d))
							f+=gapFilteringEffTable.get(m-d)*(1/aaSet.size());//.getProbability(aa);
					}
				}
				if(f>0){
					gapFilteringEffTable.put(m, f);
				}
			}
			//System.out.println(maxGapSizeLimit);
		}
	}
	
	public int getNumGap() {return numGap;}
	public int getMaxGapSize() {return maxGapSize;}
	public int getLength() {return length;}
	
	public GappedTag(AminoAcidSet aaSet){
		initAAMassTable(aaSet);
	}
	
	public GappedTag(GappedTag gt, AminoAcidSet aaSet){
		initAAMassTable(aaSet);
		length = gt.length;
		
		aaResidues = new char[length];
		gaps = new int[length];
		isAminoAcid = new BitSet(length);
		
		
		for(int i=0; i<length; i++){
	    	gaps[i] = gt.gaps[i];
	    	aaResidues[i] = gt.aaResidues[i];
		}
		nTermMass = gt.nTermMass;
		cTermMass = gt.cTermMass;
		isAminoAcid.or(gt.isAminoAcid);
		numGap = gt.numGap;
		maxGapSize = gt.maxGapSize;
		filteringEff = gt.filteringEff;
	}
	
	public void increaseGapSizeBy(int inc, int gapindex){
		if(isAminoAcid(gapindex)) return;
		if(maxGapSize == gaps[gapindex]) maxGapSize+= inc;
		
		filteringEff /= gapFilteringEffTable.get(gaps[gapindex]);
		gaps[gapindex] += inc;
		filteringEff *= gapFilteringEffTable.get(gaps[gapindex]);
		//gapFilteringEffTable.get(gaps[i]);
	}
	
	static public ArrayList<GappedTag> getGappedTags(int[] gp, int startPosition, int length, AminoAcidSet aaSet){

		//for(int i=0; i<gp.length;i++) System.out.print(" " + gp[i]);
		assert(gp.length-length+1 > startPosition);
		
		GappedTag gt = new GappedTag(aaSet);
		
		gt.length = length;
		gt.aaResidues = new char[length];
		gt.gaps = new int[length];
		gt.isAminoAcid = new BitSet(length);
		gt.nTermMass = gp[0] - gp[startPosition];
		gt.cTermMass = gp.length-length == startPosition? 0 : gp[startPosition + length];

		for(int i=0;i<length;i++){
			gt.gaps[i] = (short)(i+startPosition+1 < gp.length ? gp[i+startPosition] - gp[i+startPosition+1] : gp[i+startPosition]);
			gt.maxGapSize = gt.maxGapSize > gt.gaps[i] ? gt.maxGapSize : gt.gaps[i];
			if(gt.isThisGapAminoAcid(gp[0]-gp[i+startPosition], gt.gaps[i])){
				gt.isAminoAcid.set(i);
				gt.aaResidues[i] = aaMassTable.get(gt.gaps[i]).get(0).getResidue();
				
				if(aaSet.getAminoAcid('p') != null)
					if(gt.aaResidues[i] == 'p') gt.aaResidues[i] = 'P';
					
				gt.filteringEff *= 1/aaSet.size();//aaSet.getProbability(aaMassTable.get(gt.gaps[i]).get(0));
			}else{
				if(gt.gaps[i] > gt.maxGapSizeLimit)
					gt.filteringEff *=gapFilteringEffTable.get((short)gt.maxGapSizeLimit);
				else gt.filteringEff *= gapFilteringEffTable.get(gt.gaps[i]);
			}
		}

		gt.numGap = length - gt.isAminoAcid.cardinality();
		
		ArrayList<GappedTag> ret = new ArrayList<GappedTag>();
		
		GappedTag gt_after = new GappedTag(gt, aaSet);
		ret.add(gt_after);
		
		for(int i=0;i<length;i++){
			if(gt.isAminoAcid(i)){
				for(int j=1;j<aaMassTable.get(gt.gaps[i]).size();j++){
					for(int k=0; k< ret.size(); k++){
						GappedTag g = ret.get(k);
						gt_after = new GappedTag(g, aaSet);
						gt_after.aaResidues[i] = aaMassTable.get(gt_after.gaps[i]).get(j).getResidueStr();
						if(aaSet.getAminoAcid('p') != null)
							if(gt_after.aaResidues[i] == 'p') gt_after.aaResidues[i] = 'P';
						
						if(!ret.contains(gt_after)){
							ret.add(gt_after); k++;
						}
						
					}
					
					
				}
				
			}
		}
		
		return ret;
	}
	
	
	public int getIndexofGap(){
		return isAminoAcid.nextClearBit(0);
	}
	
	public int getNTermMass(){ return nTermMass;}
	public int getCTermMass(){ return cTermMass;}
	public float getfilteringEff() { return filteringEff;}
	
	public boolean isMatchedto(int[] gp){
		for(int i=0; i<gp.length - length + 1; i++){
			int prefixMass = gp[0] - gp[i];
			if(prefixMass < nTermMass) continue;
			if(prefixMass == nTermMass &&
				getGappedTags(gp, i, length, aaSet).contains(this)) return true;
			//	if(this.equals(new GappedTag(gp, i, length, aaSet))) return true;
			
			if(prefixMass > nTermMass) break;
		}
		return false;
	}
	
	public boolean isMatchedto(Peptide p){
		int prefixMass=0, suffixMass=0, i=0;
		while(prefixMass < nTermMass) prefixMass+=p.get(i++).getNominalMass();
		
		if(prefixMass == nTermMass){
			int k=i;
			for(int j=0; j<length;j++){
				if(this.isAminoAcid(j)){
					if(p.get(k++).getNominalMass() == (int)gaps[j]) continue;
					else return false;
				}else{
					int mass = 0;
					while(mass < (int)gaps[j]) mass += p.get(k++).getNominalMass();
					if(mass != (int) gaps[j]) return false;
				}
			}
			for(int j=k;j<p.size();j++) suffixMass+= p.get(j).getNominalMass();
			if(suffixMass == cTermMass) return true;
		}
		
		
		return false;///
	}
	
	
	public int getGapMass(int n){
		assert(n<length);
		return (int)gaps[n];
	}
	
	public char getResidue(int n){
		assert(n<length && isAminoAcid(n));
		return aaResidues[n];
	}
	
	
	
	public boolean equals(Object other) {
	    // Not strictly necessary, but often a good optimization
	    if (this == other)
	      return true;
	    if (!(other instanceof GappedTag))
	      return false;
	    GappedTag otherA = (GappedTag) other;
	    if(otherA.length != length) return false;
	    
	    for(int i=0; i<length; i++){
	    	
	    	if(isAminoAcid(i)){
	    		if(aaResidues[i] != otherA.aaResidues[i]) return false;
	    	}else{
	    		if(gaps[i] != otherA.gaps[i]) return false;
	    	}
	    }
	    
	    return 
	      (nTermMass == otherA.nTermMass) && 
	      (cTermMass == otherA.cTermMass) && 
	      (isAminoAcid.equals(otherA.isAminoAcid));     
	  }
	 
	 public boolean isAminoAcid(int n){
		 return isAminoAcid.get(n);
	 }
	 
	 public String toString(){
		 String ret = "[" + nTermMass + "]";
		 for(int i=0; i<length; i++){
			 if(isAminoAcid(i)) ret+=aaResidues[i];
			 else ret += "[" + gaps[i] + "]";
		 }
		 ret += "[" + cTermMass + "]";
		 return ret;
	 }
	 
	 
	 private boolean isThisGapAminoAcid(int prefixMass, int gap){
			if(!aaMassTable.containsKey(gap)) return false;
			boolean isPossibleMassesInHubSet = true;
			switch(gap){
			case 113:
				if(aaSet.getAminoAcid('p') != null) isPossibleMassesInHubSet = false;
				break;
			case 114: //GG N
				if(!intHubSet.contains((short)(prefixMass-57))) isPossibleMassesInHubSet = false;
				break;
			case 128: // AG GA Q K
				if(!intHubSet.contains((short)(prefixMass-57))) isPossibleMassesInHubSet = false;
				if(!intHubSet.contains((short)(prefixMass-71))) isPossibleMassesInHubSet = false;
				break;
			case 156: // VG GV R
				if(!intHubSet.contains((short)(prefixMass-57))) isPossibleMassesInHubSet = false;
				if(!intHubSet.contains((short)(prefixMass-99))) isPossibleMassesInHubSet = false;
				break;
			case 186: // EG DA VS SV AD GE W
				if(!intHubSet.contains((short)(prefixMass-57))) isPossibleMassesInHubSet = false;
				if(!intHubSet.contains((short)(prefixMass-71))) isPossibleMassesInHubSet = false;
				if(!intHubSet.contains((short)(prefixMass-87))) isPossibleMassesInHubSet = false;
				if(!intHubSet.contains((short)(prefixMass-99))) isPossibleMassesInHubSet = false;
				if(!intHubSet.contains((short)(prefixMass-115))) isPossibleMassesInHubSet = false;
				if(!intHubSet.contains((short)(prefixMass-129))) isPossibleMassesInHubSet = false;
				break;
			}
			return isPossibleMassesInHubSet;
		}
}
