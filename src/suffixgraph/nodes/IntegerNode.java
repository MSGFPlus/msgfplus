package suffixgraph.nodes;

import java.util.ArrayList;

import msutil.AminoAcid;


/**
 * The definition of a node in the trie. No duplicated edges are allowed in
 * a node. This object is used from the graph construction
 * @author jung
 *
 */
public class IntegerNode extends Node {

  
  // inner class edge comparator
  static EdgeRuler integerEdgeRuler = new EdgeRuler() {
    @Override
    public int compareEdges(int e1, int e2) {
      return e1 - e2;
    }

    @Override
    public int fromAa(AminoAcid aa) {
      return aa.getNominalMass();
    }

    @Override
    public float toFloat(int e) {
      return e;
    }

    @Override
    public Node nodeFactory() {
      return new IntegerNode();
    }

    @Override
    public Node nodeFactory(long[] edges, int edgeCount, int[] positions) {
      return new IntegerNode(edges, edgeCount, positions);
    }
  };
  
  
  /**
   * Default constructor. Only factory can instantiate this class
   */
  public IntegerNode() {}
  
  
  /**
   * This node can only be created through the factory
   * @param edges the edges as a long array.
   * @param edgeCount the number of edges filled by the array.
   * @param positions the start positions representing this node.
   */
  private IntegerNode(long[] edges, int edgeCount, int[] positions) {
    this.edgeCount = edgeCount;
    this.edges = edges;
    this.positions = positions;
  }
  
  @Override
  public void addEdge(int key, int nodeId) {
    int insertIndex = integerBinarySearch(this.edges, key, 0, getEdgeCount());
    addEdge(key, nodeId, (insertIndex<0) ? -insertIndex-1 : insertIndex);
  }
  
 
  @Override
  public int getEdge(int key) {
    return integerBinarySearch(this.edges, key, 0, getEdgeCount());
  } 
  
  
  @Override
  public ArrayList<Integer> getEdges(int keyMass) {
    int matchIndex = integerBinarySearch(this.edges, keyMass, 0, getEdgeCount());
    
    // navigate left and right for the matches
    if (matchIndex<0)   return new ArrayList<Integer>();
    
    int leftIndex = matchIndex;
    while (matchIndex>=0 && ((int)this.edges[matchIndex])==keyMass) {
      leftIndex = matchIndex;
      matchIndex--;
    }

    ArrayList<Integer> retValues = new ArrayList<Integer>();
    while (leftIndex<getEdgeCount() && ((int)this.edges[leftIndex])==keyMass) { 
      retValues.add((int)this.edges[leftIndex]);
      leftIndex++;
    }
    return retValues;
  }
  
  
  @Override
  public String toString() {
    String retString = "Total number of for modifiable edges " + getEdgeCount() + '\n';
    for (int i = 0; i < getEdgeCount(); i++) {
      long edge = edges[i];
      retString += "label: " + (int)edge + " (Integer Mass)";
      retString += " to " + (edge>>>32) + "\n";
    }
    return retString;
  }
    
  
  @Override
  public EdgeRuler getEdgeRuler() {
    return IntegerNode.integerEdgeRuler;
  }
  
  
/***** HELPER METHODS *****/   
  /**
   * Binary search method tailored for comparing edges encoding both the
   * composition (label) and the nodeId.
   * @param items the array to operate on.
   * @param mass the composition of the key to search for.
   * @param start the start index of the items array when doing the search (inclusive).
   * @param end the end index of the items array when doing the search (exclusive).
   * @return the index of the match or -i-1, where i is the insertion point.
   */
  static int integerBinarySearch(long[] items, int mass, int start, int end) {
    
    // base case
    if (end==start) {
      return -start-1;  
    }
    
    // base case
    if (end-start == 1) {
      int targetMass = (int)items[start];
      if (targetMass<mass) {
        // insert exactly at the end position
        return -end-1;
      }
      else if(mass<targetMass) {
        // insert at the start position
        return -start-1;
      }
      else {
        // they are equal
        return start;
      }
    }
    
    // recurse
    int midIndex = (end+start) / 2;    // check whether to recurse left or right
    int midMass = (int)items[midIndex];
    if (mass<midMass) {
      // recurse on the right
      return integerBinarySearch(items, mass, start, midIndex);
    }
    else if(midMass<mass) {
      // recurse on the left
      return integerBinarySearch(items, mass, midIndex+1, end);
    }
    else {
      // they middle item is equal
      return midIndex;
    }
  }
}