package suffixtree.results;

import suffixtree.Constants;


/**
 * This class parses a line from the mutated result file and creates an object
 * out of it
 * @author jung
 *
 */
public class MutatedResult implements Comparable<MutatedResult> {
  
  private long start;
  //private long end;
  private long position;
  private char original;
  private char mutation;
  private float massError;
  private float prob;
  
  
  //private boolean valid;
  private String protein;
  private String peptide;
  private String codon;
  //private String filename;
  private String line;

  
  
  
  /**
   * Constructor taking in a line from the result file
   * @param line
   */
  public MutatedResult(String line) {
    
    this.line = line;    
    
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
     * 12. Mutation position (relative)
     * 13. Original Amino Acid
     * 14. Mutated Amino Acid
     * 15. Detailed Identification
     */
    String[] tokens = line.split("\t");
    
    //System.out.println(line);
    
    // parse the fields
    //String[] filenameTokens = tokens[0].split("/");
    //this.filename = filenameTokens[filenameTokens.length-1];
    //this.fullpath = tokens[0];
    
    //this.scanNum = Integer.parseInt(tokens[1]);
    
    this.peptide = tokens[5].substring(2, tokens[5].length()-2);
    
    this.prob = Float.parseFloat(tokens[6]);
    
    this.protein = tokens[7];

    this.start = Long.parseLong(tokens[8]);
    //this.end = Long.parseLong(tokens[9]);
    
    this.massError = Float.parseFloat(tokens[10]);
        
    this.position = Long.parseLong(tokens[11]);

    this.original = tokens[12].charAt(0);
    
    this. mutation = tokens[13].charAt(0);
    
    this.codon = "---";
    
  }
  
  @Override
  public String toString() {
    return this.line;
  }

  @Override
  public int compareTo(MutatedResult o) {
    // reverse order sort. smaller first
    if (this.start > o.start) return 1;
    if (o.start > this.start) return -1;
    return 0;
  }
  
  public float getProb() { return this.prob; }
  
  public long getMutationPosition() { return this.position; }
  
  public long getStartPosition() { return this.start; }
  
  public float getMassError() { return this.massError; }
  
  public String getProtein() { return this.protein; }
  
  public String getCodon() { return this.codon; }
  
  public char getMutationAA() { return this.mutation; }
  
  public char getOriginalAA() { return this.original; }
  
  public String getPeptide() { return this.peptide; }

  // Protein annotation specific coordinates
  public int getOffset() { return Integer.parseInt(this.protein.split("\\s")[3]); }
  public int getShift() { return Integer.parseInt(this.protein.split("\\s")[1]); }
  public boolean isReversed() { return Integer.parseInt(this.protein.split("\\s")[2])==1; }
  public boolean isInsertion() { return this.original==suffixtree.Constants.EMPTY_AA; }
  public String getProteinName() { return this.protein.split("\\s")[0]; }
  
  public float delta() {
    float oMass = 0.0f, mMass = 0.0f;
    if (original!=suffixtree.Constants.EMPTY_AA) oMass = Constants.AA.getAminoAcid(original).getMass();
    if (mutation!=suffixtree.Constants.EMPTY_AA) mMass = Constants.AA.getAminoAcid(mutation).getMass();
    return this.massError + mMass - oMass;
  }
  
  
  public void setCodon(String codon) { this.codon = codon; }
}
