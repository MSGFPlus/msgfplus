package suffixtree.nodes;

import java.util.ArrayList;
import java.util.Arrays;

import suffixtree.edges.DirectedMassEdge;

/**
 * The more generalized data structure that allows keeping of longest prefix 
 * matches
 * @author jung
 *
 */
public class ComplexInternalNode extends InternalNode {

  private long[] prefixPositions; // record the ending position in the db
  private int positionIndex;
  private ComplexInternalNode parentNode;
  
  /**
   * Default constructor.
   */
  public ComplexInternalNode() {
    super();
    this.prefixPositions = new long[1];
    this.positionIndex = 0;
  }
  
  public ComplexInternalNode(DirectedMassEdge e) {
    super(e);
    this.prefixPositions = new long[1];
    this.positionIndex = 0;
  }
  
  /**
   * constructor taking an ending position that signals that this node
   * represents a query's end.
   * @param position
   */
  public ComplexInternalNode(int position) {
    super(position);
    this.prefixPositions = new long[1];
    this.positionIndex = 0;
  }
  
  /**
   * Store the coordinates of a prefix match.
   * @param start the start position in the database. Only 48 bits are significant.
   * @param extend the number of positions to extend from the start. Only 16 bits are significant.
   */
  public void addPartialMatch(long start, int extend) {
    // merge the extend value into the upper 2 bytes of the long value
    start |= ((long)extend)<<POSITION_SHIFT;
    
    /*
    for (long pos : prefixPositions) {
      if (pos==start) {
        System.out.println("Adding repeated position");
        System.exit(-8);
        return;      // this is a repeated position
      }
    }
    */
    
    // expand the array and put the new position at the end
    if (this.positionIndex>=prefixPositions.length) {
      // double the capacity
      this.prefixPositions = Arrays.copyOf(this.prefixPositions, this.positionIndex*2);
    }
    
    // put it at the index
    this.prefixPositions[this.positionIndex++] = start;
  } 
  
  
  private static final long PREFIX_MASK = 0x0000FFFFFFFFFFFFL;
  private static final int POSITION_SHIFT = 48;
  
  public void setParentNode(ComplexInternalNode n) {
    this.parentNode = n;
  }
  
  public ComplexInternalNode getParentNode() {
    return this.parentNode;
  }
  
  public int getPrefixMatchCount() {
    return this.positionIndex;
  }

  /**
   * Retrieve all partial matches from this node on, but excluding the coordinates
   * of the current node
   * @param results the coordinates list that where the results will be stored
   */
  public void retrieveAllCoordinates(ArrayList<Long> results) {
    for (int edgeIndex=0; edgeIndex<this.getDegree(); edgeIndex++) {
      ComplexInternalNode cin = (ComplexInternalNode)this.getEdgeAt(edgeIndex).getSink();
      for (int matchIndex=0; matchIndex<cin.getPrefixMatchCount(); matchIndex++) {
        results.add(cin.getCoordinatesAtIndex(matchIndex));
        cin.retrieveAllCoordinates(results);
      }
    }
  }
  
  public long getPrefixStartAtIndex(int index) {
    return decodeStartPosition(this.prefixPositions[index]);
  }
  
  public int getPrefixExtendAtIndex(int index) {
    return decodeExtension(this.prefixPositions[index]);
  }
  
  public long getPrefixEndAtIndex(int index) {
    return decodeEndPosition(this.prefixPositions[index]);
  }
  
  public long getCoordinatesAtIndex(int index) {
    return this.prefixPositions[index];
  }
  
///////// Static methods to decode the coordinates  
  public static long encodePositions(long start, long end) {
    long extension = end - start;
    return (extension<<POSITION_SHIFT) | start;
  }
  
  public static long decodeStartPosition(long number) {
    return number&PREFIX_MASK;
  }
  
  public static long decodeEndPosition(long number) {
    return decodeStartPosition(number)+decodeExtension(number);
  }
  
  public static int decodeExtension(long number) {
    return (int)(number>>>POSITION_SHIFT);
  }
  
}
