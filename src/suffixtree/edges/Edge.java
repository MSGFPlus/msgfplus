package suffixtree.edges;

import suffixtree.nodes.Node;


/**
 * The Edge of a SuffixTree. Edges are assumed to be byte labeled edges (array
 * of bytes).
 * @author jung
 *
 */
public abstract class Edge implements Comparable<Edge> {

  // the Node sink of this Edge
  private Node sink;
  
  /**
   * Retrieve the first label of this edge. 
   * @return the first byte of the label of this edge
   */
  public int getLabel() {
    return getLabelAt(0);
  }
  
  /**
   * Retrieve the label of the edge at a certain position in case that this is
   * a compressed edge.
   * @param offset the offset to consider.
   * @return the byte at the given position
   */
  public abstract int getLabelAt(int offset);
    
  
  /**
   * Gets the longest common prefix between two edges.
   * @param o the other edge object.
   * @return the longest common prefix between the objects.
   */
  public int getLongestCommonPrefixSize (Edge o) {
    int minSize = Math.min(this.size(), o.size());
    for (int i=0; i<minSize; i++) {
      if (this.getLabelAt(i) != o.getLabelAt(i)) return i;
    }
    return minSize;
  }
  
  public Edge getClone() {
    System.err.println("This method needs to be overriden to work!");
    System.exit(-1);
    return null;
  }
  
  /**
   * The length is the number of bytes of the edge.
   * @return
   */
  public abstract int length();
  
  /**
   * The size of this edge. This is the number of elements in this Edge. For 
   * example, a compressed edge has many elements in the edge, so the size is
   * the number of elements (mini edges) represented by the edge. However, for
   * a composite edge, the size is always one because the edge is being treated
   * as a single edge.
   * @return The number of bytes of the label of this edge.
   */
  public abstract int size();

  /**
   * The mass of this edge. Optional method.
   * @return
   */
  public abstract int mass();
  
  /**
   * Split the edge into 2 edges, adding an internal node to the split. The
   * original node is shortened and the other half of the edge is returned.
   * @param offset the offset index of the position to split
   * @return the second half of the split.
   */
  public abstract Edge split(int offset);
  
  
  /**
   * Getter method
   * @return the sink that this edge points to.
   */
  public Node getSink() { return sink; }
  
  
  /**
   * Setter method
   * @param sink the sink to use for this edge.
   */
  public void setSink(Node sink) { this.sink = sink; }
  
  
  /**
   * Getter method that returns the end index of an compressed edge. This 
   * operation is optional as explicit edges do not have this information. 
   * @return the end index of the FastaSequence of this edge.
   */
  public abstract int getEnd();
  
  /**
   * Getter method that returns the start index of an compressed edge. This
   * operation is optional as explicit edges do not have this information.
   * @return the start index of the FastaSequence of this edge.
   */
  public abstract int getStart();
  
  @Override
  public int compareTo(Edge o) {
    // Lexicographical ordering
    int minLength = Math.min(this.size(), o.size());
    for (int i=0; i<minLength; i++) {
      if(this.getLabelAt(i) < o.getLabelAt(i)) return -1;
      if(this.getLabelAt(i) > o.getLabelAt(i)) return 1;
    }
    
    if (this.size() < o.size()) return -1;
    if (this.size() > o.size()) return 1;
    return 0;
  }
}
