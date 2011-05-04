package suffixtree.matches;

import java.util.ArrayList;

import msutil.Peptide;

import sequences.MassSequence;
import suffixtree.Modification;

/**
 * This class represents a match to the database with 1 mutation
 * @author jung
 *
 */
public class ModMatchObject extends MatchObject {
  private int mass, modStart, modEnd;     
  
  /**
   * Convert this object into a line representation
   * @return
   */
  public String toText() {
    return String.format("%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s", 
        this.getStart(), 
        this.getEnd(), 
        this.mass, 
        this.modStart, 
        this.modEnd, 
        getQueryIndex(),
        this.getMatch(),
        this.getProtein());
  }
  
  /**
   * Constructor that builds a match object from a string printed with toText()
   * method
   * @param db the database
   * @param queries the array of queries containing the query for this match
   * @param line the string produced by toText()
   */
  public ModMatchObject(MassSequence db, ArrayList<ArrayList<Integer>> queries, String line) {
    String[] tokens = line.split("\t");
    this.setStart(Integer.parseInt(tokens[0]));
    this.setEnd(Integer.parseInt(tokens[1]));
    this.mass = Integer.parseInt(tokens[2]);
    this.modStart = Integer.parseInt(tokens[3]);
    this.modEnd = Integer.parseInt(tokens[4]);
    int queryIndex = Integer.parseInt(tokens[5]);

    setQuery(queries.get(queryIndex));
    setQueryIndex(queryIndex);
    
    this.setMatch(tokens[6]);
    this.setProtein(tokens[7]);
    
    // check correctness
    /*
    Peptide p = new Peptide(this.getMatchAsString(), suffixtree.Constants.AA);
    if (!p.isCorrect(getQuery())) {
      System.out.printf("%d-%s-%d %s is incorrect\n", this.getStart(), this.getMatchAsString(), this.getEnd(), getQuery().toString());
      System.exit(-9);
    }
    else {
      //System.out.printf("%s %s is correct\n", this.getMatchAsString(), query.toString());
    }
    */
  }
  
  
  /**
   * Default Constructor 
   * @param start
   * @param end
   * @param m
   * @param db
   * @param query
   * @param queryIndex
   */
  public ModMatchObject(long start, long end, Modification m, MassSequence db, ArrayList<Integer> query, int queryIndex) {
    
    long offset = db.getStartPosition(start);
    this.setStart((int)(start-offset));
    this.setEnd((int)(end-offset));

    char leftAA = '*', rightAA = '*';
    if (db.hasMass(start-1)) leftAA = db.getCharAt(start-1);
    if (db.hasMass(end+1)) rightAA = db.getCharAt(end+1);
    
    this.setMatch(leftAA + "." + db.getSubsequence(start, end) + "." + rightAA);
    this.setProtein(db.getAnnotation(start));
    
    this.mass = m.getMass();
    this.modStart = (int)(m.getStart()-start);
    this.modEnd = (int)(m.getEnd()-start);

    this.setQuery(query);
    this.setQueryIndex(queryIndex);
    
    // check correctness
    /*
    Peptide p = new Peptide(this.getMatchAsString(), suffixtree.Constants.AA);
    if (!p.isCorrect(query)) {
      System.out.printf("%d-%s-%d %s is incorrect\n", this.getStart(), this.getMatchAsString(), this.getEnd(), query.toString());
      System.exit(-9);
    }
    else {
      //System.out.printf("%s %s is correct\n", this.getMatchAsString(), query.toString());
    }
    */
  }
  
 
  
  @Override
  public int hashCode() {
    return (int)this.getStart();
  }
  
  @Override
  public boolean equals(Object o) {
    ModMatchObject other = (ModMatchObject)o;
    return this.getStart()==other.getStart() && 
           this.getEnd()==other.getEnd() && 
           this.mass==other.mass && 
           this.modStart==other.modStart;
  }
  

  /**
   * Create a peptide out of this match
   * @return the peptide object
   */
  public Peptide getPeptide() {
    return new Peptide(getMatchAsString());
  }
  
  @Override   // this is a more optimal implementation
  public String getMatchAsString() {
    StringBuffer sb =  new StringBuffer(super.getMatchAsString());
    sb.insert(this.modStart+1, String.format("%+d", this.mass));
    //System.out.println("Peptide string " + sb.toString());
    return sb.toString();
  }
  
  @Override
  public String getMatchAsStringWithFlankingAminoAcids() {
    StringBuffer sb =  new StringBuffer(super.getMatch());
    sb.insert(this.modStart+3, String.format("%+d", this.mass));
    //System.out.println("Peptide string " + sb.toString());
    return sb.toString();
  }
  
  
  public ArrayList<String> getAllMatches() {
    String unModSeq = super.getMatchAsString();
  	ArrayList<String> matches = new ArrayList<String>();
  	for (int i=this.modStart; i<this.modEnd; i++) {
  		StringBuffer sb =  new StringBuffer(unModSeq);
      sb.insert(i+1, String.format("%+d", this.mass));
      matches.add(sb.toString());
  	}
  	return matches;
  }
  
  @Override
  public Peptide getUnmodifiedPeptide() {
    return new Peptide(super.getMatchAsString());
  }
  
  
  public static String getSummaryHeader() {
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


  @Override
  public String getSummaryLine(String filename, 
                               int scanNum, 
                               String actMethod,
                               float pm, 
                               int charge, 
                               float offset) {
    
    StringBuffer sb = new StringBuffer();

    // 1. Filename
    sb.append(filename);
    sb.append("\t");
    
    // 2. Scan number
    sb.append(scanNum);
    sb.append("\t");
    
    // 3. Activation method
    sb.append(actMethod);
    sb.append("\t");
    
    // 4. Precursor mass
    sb.append(String.format("%.2f\t", pm));
    
    // 5. Charge
    sb.append(charge);
    sb.append("\t");
    
    // 6. The match with the modification
    sb.append(getMatchAsStringWithFlankingAminoAcids());
    sb.append("\t");
    
    // 7. Probability
    sb.append(String.format("%.1e", this.getProb()));
    sb.append("\t");
    
    // 8. The protein name
    sb.append(getProtein());
    sb.append("\t");
    
    // 9. Start
    sb.append(getStart());
    sb.append("\t");
    
    // 10. End
    sb.append(getEnd());
    sb.append("\t");
    
    // 11. Offset
    sb.append(String.format("%.3f", offset));
    
    return sb.toString();
  }
  
  
}
