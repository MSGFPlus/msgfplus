package edu.ucsd.msjava.suffixarray;

import java.util.ArrayList;
import java.util.HashMap;



/**
 * This class represents a set of matches specified as a set of (start, end) 
 * list of objects.
 * @author jung
 *
 */
public class MatchSet {
  

  
/***** HELPING INNER CLASSES *****/
  /**
   * Inner class encoding for a specified match in the set.
   * @author jung
   *
   */
  private class Match {
    private int start, end;
    
    public Match(int start, int end) {
      this.start = start;
      this.end = end;
    }
  }
  

  
/***** MEMBERS HERE *****/
  private ArrayList<Match> items; 

  
  
/***** CLASS DEFINITIONS HERE *****/
  /**
   * Default constructor.
   * @param start the index of the starting position.
   * @param end the index of the ending position.
   */
  public MatchSet() {
    this.items = new ArrayList<Match>();
  }

  
  /**
   * Add a match item to this object.
   * @param start The starting position of the match in the sequence (close interval).
   * @param end The ending position of the match (open interval).
   */
  public void add(int start, int end) {
    this.items.add(new Match(start, end));
  }
  
  
  /**
   * The number of items in this set.
   * @return the number of items in this MatchSet.
   */
  public int getSize() {
    return items.size();
  }
  
  
  /**
   * Get the starting position for the position-ith item in this set.
   * @param position
   * @return
   */
  public int getStart(int position) {
    return items.get(position).start;
  }
  
  public int getEnd(int position) {
    return items.get(position).end;
  }
  
  
  /**
   * O(n+m) intersection algorithm.
   * @param other
   * @return
   */
  public MatchSet intersect(MatchSet other) {
    // the end indexes of this object
    HashMap<Integer, Integer> ends = new HashMap<Integer, Integer>();
    MatchSet result = new MatchSet();
    
    for(Match m : this.items) {
      ends.put(m.end, m.start);
    }
    
    for(Match m : other.items) {
      Integer start = ends.get(m.start);
      if(start != null) {
        // there is a match
        result.add(start, m.end);
      }
    }
    return result;
  }
  
}
