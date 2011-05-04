package suffixtree.matches;

import java.util.ArrayList;

import msutil.Peptide;
import sequences.MassSequence;


/**
 * This is the object representation of a match for a gapped query against a
 * amino acid fasta database.
 * @author jung
 *
 */
public class ExactMatchObject extends MatchObject {

  
  /**
   * Constructor from a line serialized by the toText() method
   * @param queries the arraylist of queries to select the query for this object
   * @param line the line entry serialized by the toText() method
   */
  public ExactMatchObject(ArrayList<ArrayList<Integer>> queries, String line) {
    
    // start end queryIndex match protein. Tab separated
    String[] tokens = line.split("\t");
    
    this.setStart(Integer.parseInt(tokens[0]));
    this.setEnd(Integer.parseInt(tokens[1]));

    int queryIndex = Integer.parseInt(tokens[2]);
    this.setQuery(queries.get(queryIndex));
    this.setQueryIndex(queryIndex);

    this.setMatch(tokens[3]);
    this.setProtein(tokens[4]);
    
    // check correctness
    /*
    if (!Peptide.isCorrect(this.getMatchAsString(), query, suffixtree.Constants.AA)) {
      System.err.printf("%s %s is incorrect\n", this.getMatchAsString(), query.toString());
      System.exit(-9);
    }*/
  }
  
  
  /**
   * Default Constructor taking the complete meta information for this match
   * @param db the database
   * @param start the absolute start position of the match
   * @param end the absolute end position of the match
   * @param query the query as an arraylist of integers
   * @param queryIndex the index of this query
   */
  public ExactMatchObject(MassSequence db, long start, long end, ArrayList<Integer> query, int queryIndex) {
    
    char leftAA = '*', rightAA = '*';
    if (db.hasMass(start-1)) leftAA = db.getCharAt(start-1);
    if (db.hasMass(end+1)) rightAA = db.getCharAt(end+1);
    long offset = db.getStartPosition(start);
    
    this.setStart((int)(start-offset));
    this.setEnd((int)(end-offset));
    this.setMatch(leftAA + "." + db.getSubsequence(start, end)+ "." + rightAA);
    this.setProtein(db.getAnnotation(start));
    this.setQuery(query);
    this.setQueryIndex(queryIndex);
    
    // check correctness
    /*
    if (!Peptide.isCorrect(this.getMatchAsString(), query, suffixtree.Constants.AA)) {
      System.err.printf("%s %s is incorrect\n", this.getMatchAsString(), query.toString());
      System.exit(-9);
    }*/
  }
  
  
  
  /**
   * Return a compressed representation of the coordinates of this match based
   * on the 48-16 bit start-extension encoding of a match.
   * @return the unique coordinate representation of this match in the db.
   */
  //public long getCoors() { return ComplexInternalNode.encodePositions(this.start, this.end); }
  
 
  /*
  public ArrayList<Integer> getMatchAsArray() {
    ArrayList<Integer> retVal = new ArrayList<Integer>();
    for (long i=start; i<end; i++) {
      retVal.add(getDb().getIntegerMass(i)); 
    }
    return retVal;
  }*/
  
  /**
   * Retrieve the original query as an array of masses
   * @return the original query as an array of masses
   */
  public ArrayList<Integer> getQueryAsArray() {
    return this.getQuery();
  }
  
  @Override
  public boolean equals(Object o) {
    ExactMatchObject other = (ExactMatchObject)o;
    return this.getStart()==other.getStart() && this.getEnd()==other.getEnd();
  }
  
  @Override
  public int hashCode() {
    return this.getStart();
  }
  
  
  @Override
  public Peptide getPeptide() {
    return new Peptide(getMatchAsString());
  }



  /**
   * Outputs the header for the text representation of this object  
   * @return
   */
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
  public Peptide getUnmodifiedPeptide() {
    return getPeptide();
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
    //sb.append(getDb().getAnnotation(getStart()));
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
  
  
  
  @Override
  public String toString() {
    // start end queryIndex match protein. Tab separated
    return String.format("%d\t%d\t%d\t%s\t%s", 
        this.getStart(), 
        this.getEnd(), 
        this.getQueryIndex(),
        this.getMatch(),
        this.getProtein());
  }
}
