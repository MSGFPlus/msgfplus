package suffixtree;


/**
 * Helper class that represents a mutation
 * @author jung
 *
 */
public class Modification {
  
  private int nextStart, mass;
  long start, end; // start is inclusive, end is exclusive
  
  public Modification(int nextStart, int mass, long start, long end) {
    this.nextStart = nextStart;
    this.mass = mass;
  	this.start = start;
    this.end = end;
  }
  
  @Override
  public String toString() {
    return String.format("%d. Next start: %d", this.mass, nextStart);
  }
  
  
  public int getMass() { return this.mass; }
  
  
  /**
   * Get the coordinates of the next starting position to match.
   * @return the position to start matching next, i.e. the ending position of 
   * this mutated match plus 1.
   */
  public int getNextStart() { return this.nextStart; }
  
  public long getStart() { return this.start; }
  
  public long getEnd() { return this.end; }
}
