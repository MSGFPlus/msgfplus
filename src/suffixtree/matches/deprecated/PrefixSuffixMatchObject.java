package suffixtree.matches.deprecated;

import java.util.ArrayList;

import msutil.Peptide;

import sequences.MassSequence;
import suffixtree.Constants;
import suffixtree.matches.MatchObject;

/**
 * This class represents a prefix suffix match.
 * @author jung
 *
 */
public class PrefixSuffixMatchObject extends MatchObject {

  private long start, end; 
  private ArrayList<Long> midEnds;     // the middle ends
  private ArrayList<Long> midStarts;   // the middle starts

  
  
  /**
   * Constructor taking the prefix and suffix coordinates, the original query and
   * the database of the match. 
   * @param start the start coordinate of the prefix
   * @param midEnd the end of the prefix
   * @param midStart the start of the suffix
   * @param end the end coordinate of the suffix
   * @param db the database which this matches to
   * @param query the original query
   * @param queryIndex the index of the query in the conglomerate of sequences
   */
  public PrefixSuffixMatchObject(long start, long midEnd, long midStart, long end, MassSequence db, ArrayList<Integer> query, int queryIndex) {
    this.start = start;
    this.end = end;
    setQuery(query);
    //setDb(db);
    setQueryIndex(queryIndex);
    this.midEnds = new ArrayList<Long>();
    this.midStarts = new ArrayList<Long>();
    this.midEnds.add(midEnd);
    this.midStarts.add(midStart);
  }

  /**
   * Add alternative breaks that also make this correct
   * @param midEnd
   * @param midStart
   */
  public void addMiddle(long midEnd, long midStart) {
    this.midStarts.add(midStart);
    this.midEnds.add(midEnd);
  }
  
  
  /*
  public ArrayList<String> getAllBreakPoints() {
  	
    ArrayList<String> ret = new ArrayList<String>();
    for (int i=0; i<this.midEnds.size(); i++) {
    	long midEnd = this.midEnds.get(i)-this.getDb().getStartPosition(this.midEnds.get(i));
    	long midStart = this.midStarts.get(i)-this.getDb().getStartPosition(this.midStarts.get(i));
    	ret.add(String.format("(%d:%d)", midEnd, midStart));
    }
    return ret;
  }
  */
 

 
  
 
  
  @Override
  public Peptide getPeptide() {
    return new Peptide(getMatchAsString(), Constants.AA);
  }
  
  public ArrayList<Peptide> getPeptides() {
    ArrayList<Peptide> results = new ArrayList<Peptide>();
    for (int index=0; index<this.midEnds.size(); index++) {
      results.add(new Peptide(getMatchAsString()));
    }
    return results;
  }
  
  @Override
  public Peptide getUnmodifiedPeptide() {
    return getPeptide();
  }
  
  @Override
  public boolean equals(Object o) {
    PrefixSuffixMatchObject other = (PrefixSuffixMatchObject)o;
    if (this.start==other.start && this.end==other.end) return true;
    return false;
  }
  
  @Override
  public int hashCode() {
    return (int)this.start;  
  }
 

  public static String getSummaryHeader() {
    // TODO Auto-generated method stub
    return null;
  }

  @Override
  public String getSummaryLine(String filename, int scanNum, String actMethod,
      float pm, int charge, float offset) {
    // TODO Auto-generated method stub
    return null;
  }
  
}
