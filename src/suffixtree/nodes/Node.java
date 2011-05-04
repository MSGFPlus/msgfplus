package suffixtree.nodes;

import java.util.ArrayList;
import java.util.HashSet;

import suffixtree.edges.Edge;

/**
 * The Node of a SuffixTree allowing a fast construction of the SuffixTree.
 * @author jung
 *
 */
public interface Node {
  
  /**
   * Replace the current positions with this set of positions.
   * @param positions the positions to set this node to
   */
  public void setPositions(int[] positions);
  
  /**
   * Edge insertion method that adds an edge at the specified index. Note that 
   * there is no error check for order conservation after this operation.
   * @param edge the edge object to add.
   * @param index the index position to insert this.
   */
  public void insert(Edge edge, int index);
  
  /**
   * Replace the edge at the given index with a new edge.
   * @param index the index of the edge to replace
   * @param e the new edge
   */
  public void setEdgeAt(int index, Edge e);
  
  /**
   * Edge insertion for this node.
   * @param edge the edge to insert.
   * @return
   */
  public void insert(Edge edge);
  
  /**
   * Add stating positions to the node.
   * @param position the position index to add.
   * @return the Node after applying the given operation.
   */
  public void addPosition(int position);
  
  /**
   * Add multiple starting positions at the same time.
   * @param positions the array of positions.
   * @return the Node after applying the given operation.
   */
  public void addPositions(ArrayList<Integer> positions);
  
  /**
   * Getter method.
   * @return the stating position encoded by this node.
   */
  public int[] getPositions();
  
  /**
   * Get the number of positions (query indexes) stored in this node
   * @return the number of indexes stored by this node
   */
  public int getPositionsCount();
  
  /**
   * Search for the given edge in this node.
   * @param edge the key edge to search for.
   * @return the index of the matching edge or -i-1 if no matches were found,
   *         where i is the index of the correct insertion position.
   */
  public int search(Edge edge);
  
  /**
   * Public interface to search for an edge using a integer key  
   * @param key the key to look for
   * @return the index if match, or -index-1 of the insertion point if not matched.
   */
  public int search(int key);
  
  /**
   * Public interface to search for an edge using a integer key. Optimized method
   * that specifies the lower coordinate to use to search
   * @param key the key to look for
   * @param lower the smallest index (inclusive) to use when searching
   * @return the index if match, or -index-1 of the insertion point if not matched.
   */
  public int search(int key, int lower);
  
  /**
   * Public interface to search for an edge using a integer key. Optimized method
   * that specifies which coordinates to use to search
   * @param key the key to look for
   * @param lower the smallest index (inclusive) to use when searching
   * @param upper the largest index (exclusive) to use when searching
   * @return the index if match, or -index-1 of the insertion point if not matched.
   */
  public int search(int key, int lower, int upper);
  
  /**
   * Get the degree of this node.
   * @return the number of outgoing edges.
   */
  public int getDegree();
  
  /**
   * Get the edge at a specific index.
   * @param index the index of the edge to retrieve.
   * @return
   */
  public Edge getEdgeAt(int index);
  
  /**
   * Method to retrieve the starting position of this node and all the starting
   * position of the descendants of this node.
   * @param results
   */
  public void getAllPositions(HashSet<Integer> results);
  
  /**
   * Set the edges of this node to the given array.
   * @param edges the array of edges to use
   * @return the Node object resulting from the operation. 
   */
  public Node setEdges(Edge[] edges);
  
  /**
   * Set the edges of this node to the given array, but only use a limited 
   * number of elements from the array.
   * @param edges the array of edges to use
   * @param degree the number of elements to use
   * @return the Node object resulting from the operation.
   */
  public Node setEdges(Edge[] edges, int degree);
  
  /**
   * Retrieve the maximum edge in this node
   * @return the maximum edge in this node or null if this is a leaf node
   */
  public Edge getMaximumEdge();
  
  /**
   * Retrieve the minimum edge in this node or null if this is a leaf node
   * @return the minimum in this node
   */
  public Edge getMinimumEdge();
}
