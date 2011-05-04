package suffixtree.matches;

import java.util.ArrayList;

import suffixtree.nodes.ComplexInternalNode;
import msutil.Peptide;


/**
 * This interface presents a match object.
 * @author jung
 *
 */
public abstract class MatchObject {

  // the pointer to the original query
  private ArrayList<Integer> query;
  
  // the index of this query out of the all the queries
  private int queryIndex;
  
  // stores the probability of this match. This needs to be set
  private float p = 1f;
  
  // the annotation string with the flanking masses and the protein from the fasta file
  private String match, protein;
  
  // the relative start and end positions
  private int start, end; 
  
  
  /**
   * Create a peptide out of this match
   * @return the peptide object
   */
  public abstract Peptide getPeptide();
  
  /**
   * Return the peptide in its unmodified form if modified
   * @return
   */
  public abstract Peptide getUnmodifiedPeptide();
  
  /**
   * Check that the resulting match is indeed agreeable with the query.
   * @return the correctness value
   */
  public boolean isCorrect() {
    return Peptide.isCorrect(getPeptide().toString(), this.query, suffixtree.Constants.AA);  
  }
  
  /**
   * This is hashcode like function that compresses the start and end coordinates
   * into a long. It uses a 48-bit for the start and the rest for the offset.
   * @return the long number uniquely encoding for the coordinates in the db.
   */
  public long getEncodedCoors() {
    return ComplexInternalNode.encodePositions(getStart(), getEnd());
  }
  
  
  /**
   * Getter method of the prob of this match
   * @return the probability of this match computed by MSGF.
   */
  public float getProb() { return this.p; }
  
  /**
   * Setter method for the probability of this match
   * @param p the probability of this match
   */
  public void setProb(float p) { this.p = p; }
  
  
  /**
   * Getter method for the query index
   * @return
   */
  public int getQueryIndex() { return this.queryIndex; }
  
  /**
   * Setter method for the query index of this MatchObject
   * @param queryIndex
   */
  public void setQueryIndex(int queryIndex) { this.queryIndex = queryIndex; }
  
  
  /**
   * Setter method for the query originating this match
   * @param query the query
   */
  public void setQuery(ArrayList<Integer> query) { this.query = query; }
  
  /**
   * Getter method for the query item generating this MatchObject
   * @return the query that originated this match.
   */
  public ArrayList<Integer> getQuery() { return this.query; }
  
  
  public void setProtein(String protein) { this.protein = protein; }
  public String getProtein() { return this.protein; }
  
  public void setMatch(String match) { this.match = match; }
  public String getMatch() { return this.match; };
  
  public void setStart(int start) { this.start = start; }
  public int getStart() { return this.start; }
  
  public void setEnd(int end) { this.end = end; }
  public int getEnd() { return this.end; }
  
  


  
  /**
   * Retrieve the original representation of the match in the database with
   * flanking amino acids. If there is no flanking letter, an asterisk is used.
   * @return the original representation of the match in the database with flanking
   * amino acids.
   */
  public String getMatchAsStringWithFlankingAminoAcids() {
    return this.match;
  }
  
  /**
   * Get the flaking amino acid at the left of this sequence
   * @return the amino acid character or '*' if not applicable
   */
  public char getLeftFlankingAA() { return this.match.charAt(0); }
  
  /**
   * Represent this match as a match with coordinates and flanking matches
   * @return the string representation described above.
   */
  public String getStringWithCoordinates() { return getStart() + "-" + getMatchAsStringWithFlankingAminoAcids() + "-" + getEnd(); }
  
  /**
   * Get the flanking amino acid at the right of this sequence
   * @return the amino acid character or '*' if not applicable
   */
  public char getRightFlankingAA() { return this.match.charAt(this.match.length()-1); }
  
  /**
   * Retrieve the original representation of the match in the database.
   * @return the original representation of the match in the database.
   */
  public String getMatchAsString() { return this.match.substring(2, this.match.length()-2); }
  
 

  /**
   * Get the representation of this object in 1 line (Annotation, Protein, Score, Prob)
   * @return the tab separated line
   */
  //public abstract String getSummaryLine();
  
  
  /**
   * Get the text representation of the match object. To make things consistent
   * all child objects should call this method to serialize themselves. There
   * are some parameters required for calling this method because that information
   * is not stored in the match object. Here is the order in which the fields
   * should be printed:
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
   * @param filename filename
   * @param scanNum scan number
   * @param actMethod activation method
   * @param pm precursor mass
   * @param charge charge
   * @param offset the experimental to theoretical mass difference
   * @return the tab separated formatted string
   */
  public abstract String getSummaryLine(String filename, int scanNum, String actMethod, float pm, int charge, float offset);
  
  
  /**
   * Converts the array of integers into a string separated by commas
   * @param a the array list of integers
   * @return the string represented by the array (separated by commas)
   */
  public static String array2string(ArrayList<Integer> a) {
    StringBuffer sb = new StringBuffer();
    for (int i : a) {
      sb.append(i);
      sb.append(",");
    }
    return sb.substring(0, sb.length()-1);
  }
}
