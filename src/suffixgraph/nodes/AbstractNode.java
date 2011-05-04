package suffixgraph.nodes;

import java.util.ArrayList;

/**
 * Superclass for all the nodes.
 * @author jung
 *
 */
public abstract class AbstractNode {
  
  long[] edges;    // out edges. Upper 4 bytes are the nodeId. Lower 4 bytes are the composition.
  int[] positions; // start positions of this end node

  
  /**
   * Getter method for the starting positions that lead to this node.
   * @return the list of starting positions in the Fasta sequence.
   */
  public int[] getPositions() {
    return (this.positions==null) ? new int[0] : this.positions;
  }
  
  
  /**
   * Retrieve the edge with the given index. Note the encoding of the edges.
   * @param index the index of the element to retrieve.
   * @return the long with left 4 bytes as the nodeId and right 4 bytes as 
   *         the composition.
   */         
  public long getEdgeAt(int index) {
    return this.edges[index];  
  }
  
  /**
   * Get the nodeId encoded in the indicated edge.
   * @param index the index of the edge of interest.
   * @return the nodeId of the edge with the given index.
   */
  public int getNodeIdAt(int index) {
    return (int)(this.edges[index]>>>32);
  }
  
  /**
   * Retrieve the edge index that matches the given key edge.
   * @param key the node with the key we are looking for.
   * @return the index of the matching edge, or -i-1 where i is the index where
   *         the match would have been if it had existed. This is consistent 
   *         with return values of the binary search function.
   */
  public abstract int getEdge(int key);

  /**
   * Get the label of the edge at a certain position.
   * @param index the index of the item to retrieve.
   * @return the edge composition at the index.
   */
  public int getKeyAt(int index) {
    return (int)(this.edges[index]);
  }
  
  /**
   * Retrieve all the edge indices that match the given key edge.
   * @param keyComp the composition to look for.
   * @return the list of matching edges or an empty list if no matches were found.
   */
  public abstract ArrayList<Integer> getEdges(int keyComp);
  
  /**
   * Getter method.
   * @return the number of edges in this node.
   */
  public abstract int getEdgeCount();
  
  /**
   * Returns an object that defines the order of the nodes and other conversion 
   * operation
   * @return a comparator that defines order.
   */
  public abstract EdgeRuler getEdgeRuler();
}
