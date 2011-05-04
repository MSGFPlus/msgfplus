package suffixtree.matches;

import java.util.ArrayList;

import msutil.Peptide;

import sequences.MassSequence;
import suffixtree.Constants;
import suffixtree.Mutation;

/**
 * This class represents a match to the database with 1 mutation
 * @author jung
 *
 */
public class MutMatchObject extends MatchObject {
  private char mutation;       // the character to insert, replace, or delete
  private boolean insertMod;   // this is the flag that indicates insertion
  private int mutPos;
  
  
  
  /**
   * Convert this object into a line representation
   * @return
   */
  @Override
  public String toString() {
    int insert = 0;
    if (this.insertMod) insert = 1;
    return String.format("%d\t%d\t%d\t%c\t%d\t%d\t%s\t%s", 
        this.getStart(), 
        this.getEnd(), 
        this.mutPos, 
        this.mutation, 
        insert, 
        getQueryIndex(),
        getMatch(),
        getProtein());
  }
  
  
  /**
   * Constructor that builds an object from the result of the toText() method
   * @param db
   * @param queries
   * @param line
   */
  public MutMatchObject(MassSequence db, ArrayList<ArrayList<Integer>> queries, String line) {
    String[] tokens = line.split("\t");
    
    this.setStart(Integer.parseInt(tokens[0]));
    this.setEnd(Integer.parseInt(tokens[1]));
    this.mutPos = Integer.parseInt(tokens[2]);
    this.mutation = tokens[3].charAt(0);
    if (tokens[4].charAt(0)=='0') this.insertMod = false;
    else this.insertMod = true;
    int queryIndex = Integer.parseInt(tokens[5]);
    
    setQuery(queries.get(queryIndex));
    setQueryIndex(queryIndex);

    setMatch(tokens[6]);
    setProtein(tokens[7]);
    
    // check correctness
    /*
    if (!Peptide.isCorrect(this.getMatchAsString(), getQuery(), Constants.AA)) {
      //System.err.println(db.getCharAt(this.start-1));
      System.err.printf("%d-%s-%d (%s) %s is incorrect\n", this.getStart(), this.getMatchAsString(), this.getEnd(), this.toString(), getQuery().toString());
      System.exit(-9);
    }
    else {
      //System.err.printf("%s (%s) %s is correct\n", this.getMatchAsString(), this.toString(), query.toString());
    }*/
  }
  
  
  public MutMatchObject(long start, 
                        long end, 
                        Mutation m,
                        MassSequence db, 
                        ArrayList<Integer> query, 
                        int queryIndex) {
    
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
    
    setQuery(query);
    setQueryIndex(queryIndex);
    
    if (m.isInsertion()) {
      this.setInsertion((int)(m.getPosition()-offset), m.getMutation());
    }
    else if (m.isDeletion()) {
      this.setDeletion((int)(m.getPosition()-offset));
    }
    else {
      this.setMutation((int)(m.getPosition()-offset), m.getMutation());
    }
    
    // check correctness
    /*
    if (!Peptide.isCorrect(this.getMatchAsString(), query, Constants.AA)) {
      //System.err.println(db.getCharAt(this.start-1));
    	System.err.printf("%d-%s-%d (%s) %s is incorrect\n", this.start, this.getMatchAsString(), this.start+this.extend, this.toString(), query.toString());
    	System.exit(-9);
    }
    else {
      //System.err.printf("%s (%s) %s is correct\n", this.getMatchAsString(), this.toString(), query.toString());
    }*/
  }
  
  @Override
  public int hashCode() {
    return this.getStart();
  }
  
  public boolean isMutation() {
    if (this.insertMod==false && this.mutation!=Constants.EMPTY_AA) return true;
    return false;
  }
  
  public boolean isDeletion() {
    if (this.insertMod==false && this.mutation==Constants.EMPTY_AA) return true;
    return false;
  }
  
  public boolean isInsertion() {
    return this.insertMod;
  }
  
  
  public void setMutation(int position, char mutation) {
    this.mutPos = position;
    this.mutation = mutation;
    this.insertMod = false;
  }
  
  public void setDeletion(int position) {
    this.mutPos = position;
    this.mutation = Constants.EMPTY_AA;
    this.insertMod = false;
  }
  
  public void setInsertion(int position, char mutation) {
    this.mutPos = position;
    this.mutation = mutation;
    this.insertMod = true;
  }

  public int getMutationPosition() {
    return this.mutPos;  
  }
  
  @Override
  public boolean equals(Object o) {
    MutMatchObject other = (MutMatchObject)o;
    return this.getStart()==other.getStart() && 
           this.getEnd()==other.getEnd() && 
           this.mutPos==other.mutPos && 
           this.mutation==other.mutation;
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
    if (isDeletion()) {
      return sb.deleteCharAt(this.mutPos-this.getStart()).toString();
    }
    if (isInsertion()) {
      return sb.insert(this.mutPos-this.getStart(), this.mutation).toString();
    }
    sb.setCharAt(this.mutPos-this.getStart(), this.mutation);
    return sb.toString();
  }
  
  @Override
  public String getMatchAsStringWithFlankingAminoAcids() {
    StringBuffer sb =  new StringBuffer(super.getMatchAsStringWithFlankingAminoAcids());
    if (isDeletion()) {
      return sb.deleteCharAt(this.mutPos+2-this.getStart()).toString();
    }
    if (isInsertion()) {
      return sb.insert(this.mutPos+2-this.getStart(), this.mutation).toString();
    }
    sb.setCharAt(this.mutPos+2-this.getStart(), this.mutation);
    return sb.toString();
  }
  
  
  @Override
  public Peptide getUnmodifiedPeptide() {
    return new Peptide(super.getMatch());
  }
  
  /**
   * This is the original character in the DB
   * @return
   */
  public char getMutationSource() {
    // the was nothing in the original database
    if (this.insertMod) return Constants.EMPTY_AA;
    return super.getMatchAsString().charAt(this.getMutationPosition()-this.getStart());
  }
  
  /**
   * This is the aa that it mutates to
   * @return
   */
  public char getMutationSink() {
   return this.mutation;
  }
  

  /**
   * The mutation representation of this match with the [Sounce->Sink]
   * @return
   */
  public String toDetailedString() {
    
    StringBuffer str = new StringBuffer();
    //str.append(String.format("%d-", this.start));
    for (int index = this.getStart(); index<this.getEnd(); index++) {
      if (index != this.getMutationPosition()) {
        str.append(String.format("%c", super.getMatchAsString().charAt(index-this.getStart())));
      }
      else {
        // there is a mod here
        if (this.insertMod) {
          // insertion 
          str.append(String.format("[%c->%c]", Constants.EMPTY_AA, this.mutation));
          str.append(super.getMatchAsString().charAt(index-this.getStart()));
        }
        else {
          // substitution or deletion
          str.append(String.format("[%c->%c]", super.getMatchAsString().charAt(index-this.getStart()), this.mutation));
          //if (this.mod==Constants.EMPTY_AA) index++;
        }
      }
    }
    
    // we can have the mutation at the end
    if (this.getMutationPosition() == this.getEnd() && this.isInsertion())
      str.append(String.format("[%c->%c]", Constants.EMPTY_AA, this.mutation));
    
    //str.append(String.format("-%d", this.end));
    
    return str.toString();
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
     * 12. Mutation position (relative)
     * 13. Original Amino Acid
     * 14. Mutated Amino Acid
     * 15. Detailed Identification
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
    sb.append("MassError\t");
    sb.append("MutationPos\t");
    sb.append("OriginalAA\t");
    sb.append("MutatedAA\t");
    sb.append("DetailedIdentification");
    
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
     * 15. Detailed peptide depicting the mutation
     */
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
    
    // 6. Mutated Peptide
    sb.append(getMatchAsStringWithFlankingAminoAcids());
    sb.append("\t");

    // 7. Probability
    sb.append(String.format("%.1e", this.getProb()));
    sb.append("\t");
    
    // 8. protein identification
    sb.append(getProtein());
    sb.append("\t");
    
    // 9. Start
    sb.append(getStart());
    sb.append("\t");
    
    // 10. End
    sb.append(getEnd());
    sb.append("\t");
    
    // 11. Mass error
    sb.append(String.format("%.3f\t", offset));
    
    // 12. Mutation position
    sb.append(getMutationPosition());
    sb.append("\t");
    
    // 13. Original aa
    sb.append(getMutationSource());
    sb.append("\t");

    // 14. Mutated aa
    sb.append(getMutationSink());
    sb.append("\t");
    
    // 15. Protein identification
    sb.append(toDetailedString());
    
    return sb.toString();
  }
  
  
}
