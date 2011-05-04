package suffixgraph.nodes;

import java.util.HashSet;


/**
 * A class that holds the information about merged nodes and the merging 
 * members
 * @author jung
 *
 */
public class MergedNodeInfo {
  private int nodeId;
  private int[] items;
  
  /**
   * Default constructor taking the merged and merging id's of the nodes
   * @param nodeId the id of the merged node
   * @param items the id's of the merging nodes
   */
  public MergedNodeInfo(int nodeId, int[] items) {
    this.nodeId = nodeId;
    this.items = items;
  }
  
  
  /**
   * Getter method
   * @return the id's of the merging node
   */
  public int[] getItems() {
    return this.items;
  }
  
  
  /**
   * Getter methods
   * @return the id of the merged node
   */
  public int getNodeId() {
    return this.nodeId;
  }
  
  
  public boolean isContained(int[] others) {
    // inefficient implementation
    HashSet<Integer> tempNodeIds = new HashSet<Integer>();
    for(int other : others) tempNodeIds.add(other);
    return isContained(tempNodeIds);
  }
  
  public boolean isContained(HashSet<Integer> others) {
    for(int myItem: items) if (!others.contains(myItem)) return false; 
    return true;  
  }
  
}
