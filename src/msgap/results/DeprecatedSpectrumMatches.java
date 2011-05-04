package msgap.results;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;

import suffixtree.Constants;
import suffixtree.matches.MatchObject;

import msgap.Parameters;
import msgap.ScoringParameter;
import msgf.GeneratingFunction;
import msgf.NominalMassFactory.NominalMass;
import msutil.Peptide;



/**
 * This class represents the denovo sequenced spectrum.
 * @author jung
 *
 */
public class DeprecatedSpectrumMatches {
  // Spectrum object

  
  // The series of gapped peptide reconstructions
  //private ArrayList<ArrayList<Integer>> sequences;
  
  // The list of matches
  private HashSet<MatchObject> matches;
  
  private ScoringParameter spar = null;
  private Parameters par = null;
  
  /**
   * Constructor - initiate variables and calculate the generating function of this spectrum
   * @param s the spectrum
   */
  public DeprecatedSpectrumMatches(ScoringParameter spar, Parameters par) {
    this.matches = new HashSet<MatchObject>();
    this.par = par;
    this.spar = spar;
  }
  

  /**
   * Add the match objects. Compute the scores as the matches are added to this
   * object. This function will also automatically prune repeating matches and
   * only keep those unique matches.
   * @param index the position to set the matches
   * @param matches the matches found in the database for the sequence at the given index
   */
  public void addMatches(ArrayList<MatchObject> matches) {

    for (MatchObject mo : matches) {
      
      if (!this.matches.contains(mo)) {
        // compute the probability and scores if doable
        if (spar.getScoredSpec()!=null) {
          Peptide pep = new Peptide(mo.getMatchAsString(), par.aaSet());
          int score = getScore(pep);
          //mo.setScore(score);
          mo.setProb(getProbability(pep, mo.getMatchAsStringWithFlankingAminoAcids(), score));
        }
        this.matches.add(mo);
      }
    }
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
  private int getScore(Peptide pep){
    if(!par.allowNonEnzymaticCleavage() && !par.enzyme().isCleaved(pep)) {
      return spar.getScoredSpec().getScoreMinThreshold();
    }
    
    if(Math.abs(pep.getParentMass() - spar.getOriginalParentMass()) > par.pmTolerance().getToleranceAsDa(spar.getOriginalParentMass())) {
      //System.out.println(scoredSpec);
      return spar.getScoredSpec().getScoreMinThreshold();
    }
      
    return spar.getGraph().getScore(pep); // modified by Sangtae
  }
  
  /**
   * Gets the probability of the input peptide when peptide score is given
   * @param pep peptide 
   * @param psmScore peptide score
   * @return the spectral probability
   */
  private float getProbability(Peptide pep, String id, int psmScore){
    if(psmScore == spar.getScoredSpec().getScoreMinThreshold()) return 1;
      
    float specProb;
    
    GeneratingFunction<NominalMass> gf = new GeneratingFunction<NominalMass>(spar.getScoredSpec(), spar.getGraph()).enzyme(par.enzyme());
    
    if(par.enzyme().isCleaved(pep))
      specProb = spar.getFinalDistWellCleaved().getSpectralProbability(psmScore);
    else
      specProb = spar.getFinalDist().getSpectralProbability(psmScore);
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
   * Get the total number of matches (MatchObject) stored in this datastructure.
   * @return the count of MatchObjects
   */
  /*
  public int getMatchesCount() {
    return this.matches.size();  
  }*/
  
  
  /**
   * Gets the header fields for the text output of instance objects of this class
   * @return the header fields as a String
   */
  public static String getHeader() {
    StringBuffer sb = new StringBuffer();
    
    sb.append("SpecFile");
    sb.append("\t");
    
    sb.append("ScanNumber");
    sb.append("\t");
    
    sb.append("ActMethod");
    sb.append("\t");
    
    sb.append("Precursor");
    sb.append("\t");
    
    sb.append("Charge");
    sb.append("\t");
    
    sb.append("Annotation");
    sb.append("\t");
    
    sb.append("Protein");
    sb.append("\t");
    
    sb.append("PeptideScore");
    sb.append("\t");
    
    sb.append("SpecProb");
    sb.append("\t");
    
    return sb.toString();
  }
  
}
