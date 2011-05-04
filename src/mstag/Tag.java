/**
 * 
 */
package mstag;

import msutil.*;


/**
 * A Tag is an incomplete annotaion of a Spectrum.
 * @author jung
 *
 */
public class Tag extends Sequence<Mass> implements Comparable<Tag> {
  
  //this is recommended for Serializable objects
  static final private long serialVersionUID = 1L;
  
  // left and right flanking masses
  private Mass leftMass;
  private Mass rightMass;
  private Peptide sequence;
  private float score;
  
  // the list of peaks supporting this tag
  private Peak[] peaks;
    

  public Tag(Mass left, Peptide sequence, Mass right, Peak[] peaks, float score) {
    
    add(left);
//    for(Mass mass : sequence)     add(mass);
    for(AminoAcid aa : sequence)
    	add(new Mass(aa.getMass()));
    add(right);
    
    this.leftMass = left;
    this.sequence = sequence;
    this.rightMass = right;
    this.peaks = peaks;
    this.score = score;
  }
  
  
  public float getLeftMass() {
    return this.leftMass.getMass();
  }
    
 
  
  
  public String getTagStr() {
    return this.sequence.toString();
  }
  
  
  /**
   * Reverse other compare, sorted by score.
   */
  public int compareTo(Tag other) {
    if(score > other.score)       return -1;
    if(other.score > score)       return 1;
    return 0;
  }

  
  public String toString() {
    String retStr = "";
    retStr += "[" + leftMass.toString() + "]";
    retStr += sequence;
    retStr += "[" + rightMass.toString() + "]";
    return retStr + "\t" + score;
  }
  
  public float getRightMass() {
    return this.rightMass.getMass();
  }
  
  
  public float[] getMasses() {
    return null;
  }
  
  public static void main(String[] args) {
    //Tag t = new Tag();
  }
}
