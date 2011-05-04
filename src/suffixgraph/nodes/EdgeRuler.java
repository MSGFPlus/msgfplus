package suffixgraph.nodes;

import msutil.AminoAcid;

/**
 * Defined the order of two edges represented as integer and the different
 * transformations possible.
 * @author jung
 *
 */
public interface EdgeRuler {
  
  /**
   * The method that defines the order.
   * @param e1 one edge.
   * @param e2 the other edge.
   * @return positive number if e1 > e2, negative if e1 < e2, 0 otherwise.
   */
  public int compareEdges(int e1, int e2);
  
  
  /**
   * Convert the edge representation from integer to float.
   * @param e
   * @return
   */
  public float toFloat(int e);
  
  
  /**
   * Convert from the amino acid to an integer representation
   * @param a
   * @return
   */
  public int fromAa(AminoAcid a);
  
  
  /**
   * Just create a node suitable for housing this type of edge
   * @return
   */
  public Node nodeFactory();
  
  
  /**
   * Create a node suitable for housing this type of edge
   * @param edges the array of edges
   * @param edgesCount the number of valid edges in the array
   * @param positions the starting positions ending in this node
   * @return
   */
  public Node nodeFactory(long[] edges, int edgeCount, int[] positions);
}
