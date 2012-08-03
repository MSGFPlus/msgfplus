package edu.ucsd.msjava.msutil;

import java.util.ArrayList;

/**
 * This class construct a peptide representation based on a string of amino 
 * acids. The strings allows for brackets that groups amino acids together
 * effectively allowing for gaps. This data structure is backed by an array 
 * of compositions as an internal data structure.
 * @author jung
 *
 */
public class GappedPeptide extends Sequence<Composition> {
  //this is recommended for Serializable objects
  static final private long serialVersionUID = 1L;

  private ArrayList<Composition> compositions;
  private String sequence;
  private int count;        // number of letters in this data structure
  
  /**
   * Construct a peptide based on the string. Brackets are allowed for
   * grouping of amino acids. For example, [AC]F[DR]G.
   * @param sequence
   */
  public GappedPeptide(String sequence) {
    
    this.sequence = sequence;
    this.compositions = new ArrayList<Composition>();
    this.count = 0;
    
    boolean inBracket = false;
    Composition current = new Composition(0);
    
    for (int index=0; index<sequence.length(); index++) {
      char thisChar = sequence.charAt(index);
      if (thisChar == '[') {
        inBracket = true;
        current = new Composition(0);
      }
      else if (thisChar == ']') {
        inBracket = false;
        this.compositions.add(current);
      } 
      else {
        this.count++;
        if (inBracket) {
          current = current.getAddition(AminoAcid.getStandardAminoAcid(thisChar).getComposition());
        }
        else {
          this.compositions.add(AminoAcid.getStandardAminoAcid(thisChar).getComposition());
        }
      }
    }
  }
  
  
  /**
   * Returns the number of letters in this gapped peptide.
   * @return number of amino acids in this peptide.
   */
  public int getCount() {
    return count;
  }
  
  
  @Override
  public int size() {
    return this.compositions.size();
  }
  
  
  @Override
  public Composition get(int index) {
    return compositions.get(index);
  }
  
  
  /**
   * Return the compositions of this object.
   * @return the compositions of this object.
   */
  public ArrayList<Composition> getCompositions() {
    return this.compositions;  
  }
  
  
  /**
   * String representation of this object.
   * @return the string representation of this object.
   */
  public String toString() {
    return this.sequence;
  }

}
