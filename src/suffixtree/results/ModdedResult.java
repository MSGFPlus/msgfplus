package suffixtree.results;

import msutil.Composition;
import msutil.Peptide;


public class ModdedResult implements Comparable<ModdedResult> {
  
  private long start;
  private long end;
  private float offset;
  private float prob;
  private String protein;
  private String peptide;
  private String filepath;
  private int scanNum;
  //private String actMethod;
  private float precursorMass;
  private int charge;
  private String line;
  
  public ModdedResult(String line) {
    
    String[] tokens = line.split("\t");
    
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
    
    //System.out.println(line);
    
    //this.line = line;
    this.filepath = tokens[0];
    this.scanNum = Integer.parseInt(tokens[1]);
    //this.actMethod = tokens[2];
    this.precursorMass = Float.parseFloat(tokens[3]);
    this.charge = Integer.parseInt(tokens[4]);
    this.peptide = tokens[5];
    this.prob = Float.parseFloat(tokens[6]);
    this.protein = tokens[7];
    this.start = Long.parseLong(tokens[8]);
    this.end = Long.parseLong(tokens[9]);
    this.offset = Float.parseFloat(tokens[10]);

    this.line = line;
  }
  
  public String getFilename() {
    String[] filenameTokens = this.filepath.split("/");
    return filenameTokens[filenameTokens.length-1];
  }

  /**
   * Get the integer offset of this annotation
   * @return
   */
  public int getIntegerOffset() {
    StringBuffer sb = new StringBuffer();
    int sign = 1;
    for (int i=0; i<this.peptide.length(); i++) {
      char c = this.peptide.charAt(i); 
      if (!Character.isDigit(c)) {
        if (sb.length()>0) {
          return sign * Integer.parseInt(sb.toString());
        }
      }
      else {
        sb.append(c);
      }
      if (c=='-') sign = -1;
    }
    return 0;
  }
  
  
  public int getCharge() { return charge; }
  
  public String getModificationPositionKey() {
    return getModPosition() + "$$" + this.protein;
  }
  
  public long getModPosition() {
    for (int i=0; i<this.peptide.length(); i++) {
      char c = this.peptide.charAt(i);
      if (c=='-' || c=='+') return start+i-2; 
    }
    return end;
  }
  
  public String getPeptide() { return this.peptide; }
  
  /**
   * Calculates the exact difference between the theoretical mass and the
   * given precursor mass
   * @return
   */
  public float getDelta() {
    // create the peptide without the modification
    StringBuffer sb = new StringBuffer();
    String pep = this.peptide.split("\\.")[1];
    for (int i=0; i<pep.length(); i++) {
      if (Character.isLetter(pep.charAt(i))) {
        sb.append(pep.charAt(i));
      }
    }
    float tMass = new Peptide(sb.toString()).getMass();
    return (float)((this.precursorMass-Composition.H)*this.charge - tMass - Composition.H2O);
  }
  
  
  public String getFilepath() { return this.filepath; }
  public float getProb() { return this.prob; }
  public int getScanNumber() { return this.scanNum; }
  
  @Override
  public String toString() {
    return line;
  }

  @Override
  public int compareTo(ModdedResult o) {
    // reverse order sort. smaller first
    if (this.start > o.start) return 1;
    if (o.start > this.start) return -1;
    return 0;
  }
  
  
  public float getMassError() { return this.offset; }
  
  

}
