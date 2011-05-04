package suffixgraph.nodes;

import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;

import msutil.AminoAcid;
import msutil.Composition;


/**
 * This is the parent class for nodes that are intermediate nodes for building
 * the graphs. Defines the interfaces for adding and modifying the edges.
 * @author jung
 *
 */
public class Node extends AbstractNode {
  
  // we only need 1 instance of this per class
  static EdgeRuler compositionEdgeRuler = new EdgeRuler() {
    @Override
    public int compareEdges(int e1, int e2) {
      return Composition.compareCompositions(e1, e2);
    }

    @Override
    public float toFloat(int e) {
      return Composition.getMonoMass(e);
    }

    @Override
    public int fromAa(AminoAcid a) {
      // TODO Auto-generated method stub
      return a.getComposition().getNumber();
    }

    @Override
    public Node nodeFactory() {
      return new Node();
    }

    @Override
    public Node nodeFactory(long[] edges, int edgeCount, int[] positions) {
      return new Node(edges, edgeCount, positions);
    }
  };
  
  int edgeCount;

  
  /**
   * Default constructor.
   */
  public Node() {
    this.edgeCount = 0;
    this.edges = new long[0];
    this.positions = null;
  }
  
  
  /**
   * Constructor taking the members of this class. Only factory can instantiate this class this way.
   * @param edges the sorted edges encoded as longs.
   * @param edgesCount the number of valid edges in the array of edges.
   * @param positions the starting positions of the paths ending here.
   */
  private Node(long[] edges, int edgesCount, int[] positions) {
    this.edgeCount = edgesCount;
    this.edges = edges;
    this.positions = positions;
  }
  

  /**
   * Add an outgoing edge to this node keeping the sorted order of the edges 
   * invariant. This is a rather expensive operation because the array needs
   * to be copied into a new array, but memory efficient.
   * @param edge the edge encoding the composition and nodeId.
   * @param insertIndex the index to be inserted to keep the sorted order.
   */
  public void addEdge(long edge, int insertIndex) {
    
    if (this.edgeCount>=this.edges.length) {
      // double the capacity because can't fit anything else
      long[] tempEdges = new long[Math.max((int)this.edges.length*2,1)];   
      System.arraycopy(this.edges, 0, tempEdges, 0, insertIndex);
      System.arraycopy(this.edges, insertIndex, tempEdges, insertIndex+1, this.edgeCount-insertIndex);
      tempEdges[insertIndex] = edge;
      this.edges = tempEdges;
    }
    else {
      System.arraycopy(this.edges, insertIndex, this.edges, insertIndex+1, this.edgeCount-insertIndex);
      this.edges[insertIndex] = edge;
    }
    this.edgeCount++;
  }
  
  
  /**
   * Add an outgoing edge to this node keeping the sorted order of the edges 
   * invariant. This is a rather expensive operation because the array needs
   * to be copied into a new array, but memory efficient.
   * @param key the key of the edge (label).
   * @param nodeId the index of sink node in the main array.
   * @param insertIndex the index to be inserted to keep the sorted order.
   */
  public void addEdge(int key, int nodeId, int insertIndex) {
    long edge = (((long)nodeId)<<32) | key;
    addEdge(edge, insertIndex);
  }
  
  
  /**
   * Add an outgoing edge to this node keeping the sorted order of the edges 
   * invariant. This is a rather expensive operation because the array needs
   * to be copied into a new array, but memory efficient.
   * @param key the key of the edge (label).
   * @param nodeId the index of sink node in the main array.
   */
  public void addEdge(int key, int nodeId) {
    int insertIndex = binarySearch(this.edges, key, Composition.getMonoMass(key), 0, this.edgeCount);
    addEdge(key, nodeId, (insertIndex<0) ? -insertIndex-1 : insertIndex);
  }
  
  
  /**
   * Add a starting position for this ending node.
   * @param pos the index of the starting position terminating in this node.
   */
  public void addPosition(int pos) {
    if (this.positions==null) {
      this.positions = new int[1];
      this.positions[0] = pos;
      return;
    }
    
    // create a temp array, copy the contents and reassign the pointer. Simulate an append.
    int[] temp = new int[this.positions.length+1];
    temp[temp.length-1] = pos;
    System.arraycopy(this.positions, 0, temp, 0, this.positions.length);
    this.positions = temp;
  }
  
  
  /**
   * Replace an edge at an index with the given nodeId and composition. 
   * @param key the new composition to use.
   * @param nodeId the new nodeId to use.
   * @param index the index of the item to replace.
   */
  public void setEdge(int key, int nodeId, int index) {
    long edge = (((long)nodeId)<<32) | key;
    this.edges[index] = edge;
  }

  
  /**
   * Setter method.    
   * @param positions the new positions to update.
   */
  public void setPositions(int[] positions) {
    this.positions = positions;
  }
  
  
  /**
   * Write the given node in binary form to the given stream.
   * @param n the node.
   * @param out the stream.
   * @throws IOException if error occurs.
   */
  public void toFile(DataOutputStream out) throws IOException {
    out.writeInt(edgeCount);
    for (int i=0; i<edgeCount; i++) {
      out.writeLong(edges[i]);
    }
    if (this.positions==null) 
      out.writeInt(0);
    else {
      out.writeInt(positions.length);
      for (int i=0; i<positions.length; i++) {
        out.writeInt(positions[i]);
      }
    }
  }
  
  
  @Override
  public int getEdge(int key) {
    return binarySearch(this.edges, key, Composition.getMonoMass(key), 0, getEdgeCount());
  } 
  
  
  @Override
  public ArrayList<Integer> getEdges(int keyComp) {
    int matchIndex = binarySearch(this.edges, keyComp, Composition.getMonoMass(keyComp), 0, getEdgeCount());
    
    // navigate left and right for the matches
    if (matchIndex<0)   return new ArrayList<Integer>();
    
    int leftIndex = matchIndex;
    while (matchIndex>=0 && ((int)this.edges[matchIndex])==keyComp) {
      leftIndex = matchIndex;
      matchIndex--;
    }

    ArrayList<Integer> retValues = new ArrayList<Integer>();
    while (leftIndex<getEdgeCount() && ((int)this.edges[leftIndex])==keyComp) { 
      retValues.add((int)this.edges[leftIndex]);
      leftIndex++;
    }
    return retValues;
  }
  
  
  @Override
  public int getEdgeCount() {
    return this.edgeCount;
  }

  
  /***** HELPER METHODS *****/   
  /**
   * Binary search method tailored for comparing edges encoding both the
   * composition (label) and the nodeId.
   * @param items the array to operate on.
   * @param keyComp the composition of the key to search for.
   * @param keyMass the mass of the key to search for to avoid recomputing this.
   * @param start the start index of the items array when doing the search (inclusive).
   * @param end the end index of the items array when doing the search (exclusive).
   * @return the index of the match or -i-1, where i is the insertion point.
   */
  static int binarySearch(long[] items, int keyComp, float keyMass, int start, int end) {
    
    // base case
    if (end==start) {
      return -start-1;  
    }
    
    // base case
    if (end-start == 1) {
      int targetComp = (int)items[start];
      float targetMass = Composition.getMonoMass(targetComp);
      if (targetMass<keyMass) {
        // insert exactly at the end position
        return -end-1;
      }
      else if(keyMass<targetMass) {
        // insert at the start position
        return -start-1;
      }
      else {
        // the float masses are exactly the same, so use the composition number to break the tie
        if (targetComp<keyComp) {       
          // insert exactly at the end position
          return -end-1;
        }
        else if(keyComp<targetComp) {
          // insert at the start position
          return -start-1;
        }   
        else {  
          // they are equal
          return start;
        }
      }
    }
    
    // recurse
    int midIndex = (end+start) / 2;    // check whether to recurse left or right
    int midComp = (int)items[midIndex];
    float midMass = Composition.getMonoMass(midComp);
    if (keyMass<midMass) {
      // recurse on the right
      return binarySearch(items, keyComp, keyMass, start, midIndex);
    }
    else if(midMass<keyMass) {
      // recurse on the left
      return binarySearch(items, keyComp, keyMass, midIndex+1, end);
    }
    else {
      // float masses are exactly the same, so use the composition number to break the tie
      if (midComp<keyComp) {       
        // recurse on right
        return binarySearch(items, keyComp, keyMass, start, midIndex);
      }
      else if(keyComp<midComp) {
        // recurse on left
        return binarySearch(items, keyComp, keyMass, midIndex+1, end);
      }   
      else {  
        // they middle item is equal
        return midIndex;
      }
    }
  }


  @Override
  public EdgeRuler getEdgeRuler() {
    return compositionEdgeRuler;
  }
}
