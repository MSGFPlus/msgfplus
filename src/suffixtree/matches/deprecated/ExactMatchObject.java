package suffixtree.matches.deprecated;

import java.util.ArrayList;

import msutil.Peptide;

import sequences.MassSequence;
import suffixtree.nodes.ComplexInternalNode;


/**
 * This is the object representation of a match for a gapped query against a
 * amino acid fasta database.
 * @author jung
 *
 */
public class ExactMatchObject extends MatchObject {
  private long start;
  private long end;
  
  
  
  /**
   * Convert this object into a line representation
   * @return
   */
  public String toText() {
    return String.format("%d\t%d\t%d", this.start, this.end, getQueryIndex());
  }
  
  
  
  /**
   * Constructor from a line serialized by the toText() method
   * @param db the database
   * @param queries the arraylist of queries to select the query for this object
   * @param line the line entry serialized by the toText() method
   */
  public ExactMatchObject(MassSequence db, ArrayList<ArrayList<Integer>> queries, String line) {
    String[] tokens = line.split("\t");
    //System.out.println("Line is " + line);
    this.start = Long.parseLong(tokens[0]);
    this.end = Long.parseLong(tokens[1]);
    int queryIndex = Integer.parseInt(tokens[2]);
    
    this.setQuery(queries.get(queryIndex));
    this.setDb(db);
    this.setQueryIndex(queryIndex);
    
    // check correctness
    /*
    if (!Peptide.isCorrect(this.getMatchAsString(), query, suffixtree.Constants.AA)) {
      System.err.printf("%s %s is incorrect\n", this.getMatchAsString(), query.toString());
      System.exit(-9);
    }*/
  }
  
  
  /**
   * Default constructor.
   * @param sequence the database
   * @param start the start position of the match
   * @param end the end position of the match
   * @param query the query in which this match was originated
   */
  public ExactMatchObject(MassSequence sequence, long start, long end, ArrayList<Integer> query, int queryIndex) {
    this.setQuery(query);
    this.setDb(sequence);
    this.setQueryIndex(queryIndex);
    this.start = start;
    this.end = end;
    
    // check correctness
    /*
    if (!Peptide.isCorrect(this.getMatchAsString(), query, suffixtree.Constants.AA)) {
      System.err.printf("%s %s is incorrect\n", this.getMatchAsString(), query.toString());
      System.exit(-9);
    }*/
  }
  
  
  /**
   * Constructor that will figure out the ending position given the FastaSequence object.
   * @param sequence the database
   * @param start the start position of the match
   * @param query the query in which this match was originated
   * @param totalMass the total mass of the query
   */
  public ExactMatchObject(MassSequence sequence, long start, ArrayList<Integer> query, int totalMass, int queryIndex) {
    this.setDb(sequence);
    this.setQueryIndex(queryIndex);
    this.start = start;
    int cumMass = 0;
    for (long i=start; i<sequence.getSize(); i++) {
      cumMass += sequence.getIntegerMass(i);
      if (cumMass >= totalMass) {
        this.end = i+1;
        break;
      }
    }
    this.setQuery(query);
  }
  
  /**
   * Return a compressed representation of the coordinates of this match based
   * on the 48-16 bit start-extension encoding of a match.
   * @return the unique coordinate representation of this match in the db.
   */
  public long getCoors() {
    return ComplexInternalNode.encodePositions(this.start, this.end);  
  }
  
  @Override
  public long getStart() { return this.start; }
  
  @Override
  public long getEnd() { return this.end; }
 
  @Override
  public String getMatchAsString() {
    return getDb().getSubsequence(this.start, this.end);
  }
  
  public ArrayList<Integer> getMatchAsArray() {
    ArrayList<Integer> retVal = new ArrayList<Integer>();
    for (long i=start; i<end; i++) {
      retVal.add(getDb().getIntegerMass(i)); 
    }
    return retVal;
  }
  
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
    return this.start==other.start && this.end==other.end;
  }
  
  @Override
  public int hashCode() {
    return (int)this.start;
  }
  
  @Override
  public String toString() {
    return this.start + ":" + this.end + ":" + getMatchAsString() + ":" + array2string(this.getQuery()) + ":" + getQueryIndex();
  }
  
  public String shortSummary() {
    return this.start + ":" + this.end + ":" + getQueryIndex();
  }
  
  @Override
  public Peptide getPeptide() {
    return new Peptide(getMatchAsString());
  }

  @Override
  public String getSummaryLine() {
    StringBuffer sb = new StringBuffer();
    
    // the match in the database
    sb.append(getMatchAsStringWithFlankingAminoAcids());
    sb.append("\t");
    
    // the protein identification
    sb.append(getDb().getAnnotation(getStart()));
    sb.append("\t");
    
    // the score
    sb.append(this.getScore());
    sb.append("\t");
    
    // the probability
    sb.append(String.format("%.1e", this.getProb()));
    sb.append("\t");
    
    // the relative start and end
    sb.append(this.getRelativeStart());
    sb.append("\t");
    
    sb.append(this.getRelativeEnd());
    sb.append("\t");
    
    return sb.toString();
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
    sb.append(getDb().getAnnotation(getStart()));
    sb.append("\t");
    
    long startOffset = this.getDb().getStartPosition(getStart());
    // 9. Start
    sb.append(getStart()-startOffset);
    sb.append("\t");
    
    // 10. End
    sb.append(getEnd()-startOffset);
    sb.append("\t");
    
    // 11. Offset
    sb.append(String.format("%.3f", offset));
    
    return sb.toString();
  }
  
   
}
