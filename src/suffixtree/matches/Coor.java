package suffixtree.matches;

public class Coor {

  private long start;
  private int extend;
  
  public Coor(long start, long end) {
    this.start = start;
    this.extend = (int)(end-start);
  }
  
  public long getStart()     { return this.start; }
  public long getEnd()       { return this.start+extend; }
}
