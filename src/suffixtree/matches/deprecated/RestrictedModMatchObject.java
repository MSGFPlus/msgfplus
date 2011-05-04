package suffixtree.matches.deprecated;

import java.util.ArrayList;

import suffixtree.matches.ExactMatchObject;

import msgap.ModifiedAAinGap;
import msutil.Peptide;

// added by kyowon
public class RestrictedModMatchObject extends ExactMatchObject{
	 private Peptide modifiedPeptide;
	 
	 /**
	   * Default constructor. Use getModificationMatchObjectsFrom(ExactMatchObject emo, AminoAcidSet aaSet) to get instants
	   * @param emo the exact match that this will be based on
	   * @param modifiedPeptide which modified peptide this will be representing
	   */
	 private RestrictedModMatchObject(ExactMatchObject emo, Peptide modifiedPeptide, boolean thisClassDoesNotWork) {
		 super(null, emo.getStart(), emo.getEnd(), emo.getQuery(), emo.getQueryIndex());
		 this.modifiedPeptide = modifiedPeptide;
	 }
	 
	  @Override
	  public String getMatchAsString() {
	    return modifiedPeptide.toString();
	  }

	  @Override
	  public boolean equals(Object o) {
		  RestrictedModMatchObject other = (RestrictedModMatchObject)o;
	    return this.getStart()==other.getStart() && this.getEnd()==other.getEnd() && this.modifiedPeptide.toString().equals(other.getPeptide().toString());
	  }
	  
	  @Override
	  public Peptide getPeptide() {
	    return modifiedPeptide;
	  }

	  
	  /*The instant(s) getter*/
	  public static ArrayList<RestrictedModMatchObject> getModificationMatchObjectsFrom(ExactMatchObject emo){
		  ArrayList<RestrictedModMatchObject> mmo = new ArrayList<RestrictedModMatchObject>();
		  for(Peptide p : ModifiedAAinGap.getModifiedPeptides(emo.getPeptide()))
			  mmo.add(new RestrictedModMatchObject(emo, p, true));
		  
		  return mmo;
	  }
	  
	 
}
