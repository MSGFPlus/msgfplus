package msgap.results;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import suffixtree.matches.ExactMatchObject;
import suffixtree.matches.MatchObject;
import suffixtree.matches.deprecated.RestrictedModMatchObject;

import msgap.ModifiedAAinGap;
import msgap.Parameters;
import msgap.ScoringParameter;
import msgf.GeneratingFunction;
import msgf.NominalMass;
import msgf.ScoredSpectrum;
import msutil.Peptide;


/**
 * This class represents the denovo sequenced spectrum.
 * @author jung
 *
 */
public class SpectrumMatches {
 
  // The list of matches
  private ArrayList<MatchObject> matches;
  private ScoringParameter spar = null;
  private Parameters par = null;
  
  
  /**
   * Constructor - initiate variables and calculate the generating function of this spectrum
   * the parent mass will be corrected if correctPM in Parameters is true 
   * @param s the spectrum
   */
  public SpectrumMatches(ScoringParameter spar, Parameters par) {
    this.matches = new ArrayList<MatchObject>();
    this.par = par;
    this.spar = spar;
  }
  
  
  /**
   * Add the match objects. Compute the scores as the matches are added to this
   * object. This function will also automatically prune repeating matches and
   * only keep those unique matches.
   * @param index the position to set the matches
   * @param matches the matches found in the database for the sequence at the given index
   * @return the best probability of this set of matches
   */
  public float addMatches(ArrayList<MatchObject> matches) {
	  
	  //this part is added by kyowon for fixed mod search. it does not change anything for 
	  // normal exact matching or other matchings used by Julio
	 if(ModifiedAAinGap.isAnyAAModified() && ((matches.get(0)) instanceof ExactMatchObject)){
		 ArrayList<MatchObject> newmatches = new ArrayList<MatchObject>();
		 
		 for (MatchObject mo : matches) {
			 for(RestrictedModMatchObject mmo : RestrictedModMatchObject.getModificationMatchObjectsFrom((ExactMatchObject) mo)){
				 newmatches.add(mmo);
			 }
		 }
		 matches = newmatches;
	 }
	 //up to here, this part
	 
    float bestProb = 2.0f;
    HashMap<String,MatchObject> seen = new HashMap<String,MatchObject>();
    for (MatchObject mo : matches) {
      String flankingStr = mo.getMatchAsStringWithFlankingAminoAcids();
      String key = mo.getStart() + flankingStr;
      if (!seen.containsKey(key)) {
        // compute the probability and scores 
        Peptide pep = new Peptide(flankingStr.substring(2, flankingStr.length()-2), par.aaSet());
        int score = getScore(pep);
        //mo.setScore(score);
        float p = getProbability(pep, flankingStr, score);
        if (p!=0) {
          mo.setProb(p);
          seen.put(key, mo);
          this.matches.add(mo);
          bestProb = Math.min(bestProb, p); // register the best Prob
        }
      }
      else {
	    	// set the prob with the previous one
	    	//mo.setScore(seen.get(key).getScore());
	    	//mo.setProb(seen.get(key).getProb());
      }
    }
    
    return bestProb;
  }
  
  /**
   * Get the scan number of the spectrum stored in this object
   * @return the scan number of the spectrum in this object
   */
  public int getScanNumber() { return spar.getScanNum(); }
  

  /**
   * Gets the score of the input peptide. Parent mass tolerance considered.
   * @param pep peptide           
   * @return the peptide score
   */
  public int getScore(Peptide pep){
	//  System.out.println(this.getScanNumber() + " " + pep );
    if(!par.allowNonEnzymaticCleavage() && !par.enzyme().isCleaved(pep)) {
      return spar.getScoredSpec().getScoreMinThreshold();
    }
    
    if(Math.abs(pep.getParentMass() - spar.getCorrectedParentMass()) > par.pmTolerance().getToleranceAsDa(spar.getCorrectedParentMass())
    		|| spar.getGraph().getPMNode().getNominalMass() != pep.getNominalMass()) {
    		
    	//System.out.println(scoredSpec);
    	
      return spar.getScoredSpec().getScoreMinThreshold();
    }
      
  //  if(!par.enzyme().isCleaved(pep))
    //    spar.getGraph().allowNonEnzymaticCleavage();
        
    return spar.getScoredSpec().getNodeScore(pep, spar.getGraph()); // modified by Sangtae
  }
  
  /**
   * Gets the probability of the input peptide when peptide score is given
   * @param pep peptide 
   * @param psmScore peptide score
   * @return the spectral probability
   */
  public float getProbability(Peptide pep, String id, int psmScore){
    if(psmScore == spar.getScoredSpec().getScoreMinThreshold()) return 1;
      
    GeneratingFunction<NominalMass> gf = new GeneratingFunction<NominalMass>(spar.getScoredSpec(), spar.getGraph()).enzyme(par.enzyme());
    float specProb;
    if(par.enzyme().isCleaved(pep))
      specProb = spar.getFinalDistWellCleaved().getSpectralProbability(psmScore)*gf.getNeighboringAAProbability(id);
    else
      specProb = spar.getFinalDist().getSpectralProbability(psmScore)*gf.getNeighboringAAProbability(id);
    
  //  if(specProb <=0) System.out.println(this.getScanNumber() + " " + id + " " + psmScore + " " + specProb);
  //  assert(specProb > 0);
   //
    
    return specProb;
  }
  
  
  /**
   * Get the best probability out of all the matches in this object
   * @return the best (smallest) probability
   */
  public float getBestProbability() {
    float bestProb = Float.MAX_VALUE;
    for (MatchObject mo : this.matches) {
      if (mo.getProb() < bestProb) bestProb = mo.getProb();
    }
    return bestProb;
  }
  
  
   /**
   * Gets the probability of the input peptide
   * @param pep peptide 
   * @param id peptide string including flanking aa's     
   * @return the spectral probability
   */
  /*
  private float getProbability(Peptide pep, String id){
    return getProbability(pep, id, getScore(pep));
  }
  */

  
  public ScoredSpectrum<NominalMass> getScoredSpec() {return spar.getScoredSpec(); }
  
  /**
   * The text representation of this object in congruence with the getHeader method.
   * @return the list of probabilities.
   */
  /*
  public void textOutput(float specProbCutOff, ArrayList<String> results) {
    for (MatchObject mo : this.matches) {
      if (mo.getProb()<=specProbCutOff) {
        results.add(mo.getSummaryLine());
      }
    }
  }
  */
  
  public Collection<MatchObject> getMatches() {
    return this.matches;
  }
  
  
  /*
  @Override
  public String toString() {
    ArrayList<String> results = new ArrayList<String>();
    textOutput(Constants.PROB_CUTOFF, results);
    StringBuffer sb = new StringBuffer();
    for (String s : results) {
      sb.append(s).append("\n");
    }
    return sb.toString();
  }
  */
  
  /**
   * Gets the header fields for the text output of instance objects of this class
   * @return the header fields as a String
   */
  public static String getHeader() {
    StringBuffer sb = new StringBuffer();
    
    /**
     * 1. Filename
     * 2. Scan number
     * 3. Activation method
     * 4. Precursor mass
     * 5. Charge
     * 6. Peptide match / annotation
     * 7. Probability
     * 8. Protein name
     * 9. Start position in the protein
     * 10. End position in the protein
     * 11. Offset (Experimental Mass - Theoretical Mass)
     * 12+ Additional fields depending on the type of match object
     */
    sb.append("Filename\t");
    sb.append("ScanNum\t");
    sb.append("ActMethod\t");
    sb.append("PrecursorMass\t");
    sb.append("Charge\t");
    sb.append("Annotation\t");
    sb.append("Probability\t");
    sb.append("Protein\t");
    sb.append("Start\t");
    sb.append("End\t");
    sb.append("MassError");
    
    return sb.toString();
  }
  
}
