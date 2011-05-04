package msgap;

import java.util.ArrayList;
import java.util.Hashtable;

import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Peptide;


public class GappedTagSet extends ArrayList<GappedTag>{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private static int length = 3;
	private static int numGap = 1;
	private static int maxGapSize = 500;
	
	public static void setTagLength(int l){length = l;}
	public static void setNumGap(int n){numGap = n;}
	public static void setMaxGapSize(int m){maxGapSize = m;}
	
	public static int getTagLength(){return length;}
	public static int getNumGap(){return numGap;}
	public static int getMaxGapSize(){return maxGapSize;}
	static AminoAcidSet aaSet = null;
	//static void initAASet(AminoAcidSet a){
	//	aaSet = a;
	//}
	
	public GappedTagSet(ArrayList<int[]> deNovo, ArrayList<Integer> intHubSet, AminoAcidSet a){
		super();
		if(aaSet == null) aaSet = a;
		GappedTag.setIntHubSet(intHubSet);
		for(int[] gp : deNovo){
			if(this.isCoveringGP(gp)) continue;
			float bestFilteringEff = 1;
			int startPositionofBestTag = -1;
			for(int i=1;i<gp.length - length ; i++){ // exclude n c terminal
				GappedTag tmpTag = GappedTag.getGappedTags(gp, i, length, aaSet).get(0);
				if(tmpTag.getNumGap() > numGap || tmpTag.getMaxGapSize() > maxGapSize) continue;
				if(bestFilteringEff >= tmpTag.getfilteringEff()){ // new tag is better
					bestFilteringEff = tmpTag.getfilteringEff();
					startPositionofBestTag = i;
				}
			}
			
			if(startPositionofBestTag < 0){ // for n c terminal
				ArrayList<Integer> ncIndex = new ArrayList<Integer>();
				ncIndex.add(0); ncIndex.add(gp.length - length);
				for(int i:ncIndex){
					GappedTag tmpTag =  GappedTag.getGappedTags(gp, i, length, aaSet).get(0);
					if(tmpTag.getNumGap() > numGap || tmpTag.getMaxGapSize() > maxGapSize) continue;
					if(bestFilteringEff > tmpTag.getfilteringEff()){ // new tag is better
						bestFilteringEff = tmpTag.getfilteringEff();
						startPositionofBestTag = i;
					}
				}
			}
			
			if(startPositionofBestTag>=0){
				for(GappedTag gt : GappedTag.getGappedTags(gp, startPositionofBestTag, length, aaSet))
					if(!this.contains(gt)) this.add(gt);
			}
		}
		
		//for(GappedTag gt : this)
			//if(gt.toString().contains("K") || gt.toString().contains("Q"))
				//System.out.println("*********************************" + gt);
	}
		
	public void addGappedTagswithPTMs(){ // now only for p from P.. modify next
		//AminoAcid modifedP = aaSet.getAminoAcid('p');
		int massDiffbetweenPandp = aaSet.getAminoAcid('P').getNominalMass() - aaSet.getAminoAcid('p').getNominalMass();
		ArrayList<GappedTag> toadd = new ArrayList<GappedTag>();
		for(GappedTag gt : this){
			int max_p_num = 0;
			if(gt.getNumGap() == 0) continue; // no mass gap
			for(String aaString : GappedGeneratingFunction.getAACombTable().get(gt.getGapMass(gt.getIndexofGap()))){
				//System.out.println(aaString);
				if(aaString.contains("p")){
					int i=0;
					for(int j=0; j<aaString.length(); j++){
						char aa = aaString.charAt(j);
						if(aa == 'p') i++;
					}
					max_p_num = max_p_num < i? i : max_p_num;
				}	
			}
			//System.out.println(max_p_num);
			//System.out.println(gt+" ori");
			for(int i=1; i<max_p_num; i++){
				GappedTag newgt = new GappedTag(gt, aaSet);
				newgt.increaseGapSizeBy(massDiffbetweenPandp * i, newgt.getIndexofGap());
				//System.out.println(newgt + "\t" + newgt.getMaxGapSize() + this.contains(newgt)) ;
				if(newgt.getMaxGapSize() <= maxGapSize && !this.contains(newgt) && !toadd.contains(newgt))
					toadd.add(newgt);
			}
		//	System.out.println(this);
		}
		
		this.addAll(toadd);
		
	}
	
	public boolean isCoveringGP(int[] gp){
		for(GappedTag tag:this) if(tag.isMatchedto(gp)) return true;
		return false;
	}	
}
