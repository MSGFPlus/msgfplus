package suffixtree;


/**
 * Helper class that represents a mutation
 * @author jung
 *
 */
public class Mutation {
  
  private long position;     // absolute start position
  private long nextStart;    // where to resume matching
  private char original;
  private char mutation;
  // this flags whether the mutation was a deletion of a boundary amino acid
  // leading to a complete edge matching (the mutation was not IN the edge)
  private boolean inEdge;    
  
  public Mutation(long position, long nextStart, char original, char mutation) {
    this.position = position;
    this.nextStart = nextStart;
    this.original = original;
    this.mutation = mutation;
    this.inEdge = true;
  }
  
  public Mutation(long position, long nextStart, char original, char mutation, boolean inEdge) {
    this(position, nextStart, original, mutation);
    this.inEdge = inEdge;
  }
  
  @Override
  public String toString() {
    return String.format("%d[%c->%c]. Next start: %d. Mutation in edge: %b.", position, original, mutation, nextStart, inEdge);
  }
  
  /**
   * Check whether this mutation was inside an edge. The only case when this 
   * can happen is when there is deletion at a boundary amino acid.
   * @return whether the mutation is inside the edge
   */
  public boolean isInEdge() { return this.inEdge; }
  
  /**
   * Get the mutation represented by this object
   * @return the mutation 
   */
  public char getMutation() { return this.mutation; }
  
  
  /**
   * Get the coordinates of the next starting position to match.
   * @return the position to start matching next, i.e. the ending position of 
   * this mutated match plus 1.
   */
  public long getNextStart() { return this.nextStart; }
  
  
  /**
   * Return the coordinate of this mutation
   * @return the position of the mutation
   */
  public long getPosition() { return this.position; }
  
  
  /**
   * Check whether this mutation is a deletion.
   * @return whether this object is an deletion.
   */
  public boolean isDeletion() { return this.mutation==Constants.EMPTY_AA; }
  
  
  /**
   * Check whether this mutation is an insertion.
   * @return whether this object is an insertion.
   */
  public boolean isInsertion() { return this.original==Constants.EMPTY_AA; }
  
}
