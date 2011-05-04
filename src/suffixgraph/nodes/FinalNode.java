package suffixgraph.nodes;

import java.io.DataInputStream;
import java.io.IOException;
import java.util.ArrayList;

import msutil.Composition;


/**
 * This method is the space efficient Node that does not allow modification of 
 * its edges. We can save at least 4 bytes (Possibly 8 bytes) with this 
 * representation. The default behavior assumes that edges are comparable by 
 * composition mass. 
 * @author jung
 *
 */
public class FinalNode extends AbstractNode {

  
  /**
   * Default constructor. The only way to create an instance is through the
   * factory method.
   */
  private FinalNode() { }
  
  
  /**
   * Read a node from the input stream.
   * @param in the input stream.
   * @return a node from the input stream advancing it.
   * @throws IOException if error occurs.
   */
  public static FinalNode finalNodeFactory(DataInputStream in) throws IOException {
    FinalNode n = new FinalNode();
    int edgeCount = in.readInt();
    n.edges = new long[edgeCount];
    for (int i=0; i<edgeCount; i++)         n.edges[i] = in.readLong();
    int positionsCnt = in.readInt();
    if (positionsCnt > 0) {
      n.positions = new int[positionsCnt];
      for (int i=0; i<positionsCnt; i++)    n.positions[i] = in.readInt();
    }
    return n;
  }
   
  
  @Override
  public int getEdge(int key) {
    return Node.binarySearch(this.edges, key, Composition.getMonoMass(key), 0, getEdgeCount());
  } 
  
  
  @Override
  public ArrayList<Integer> getEdges(int key) {
    int matchIndex = Node.binarySearch(this.edges, key, Composition.getMonoMass(key), 0, getEdgeCount());
    
    // navigate left and right for the matches
    if (matchIndex<0)   return new ArrayList<Integer>();
    
    int leftIndex = matchIndex;
    while (matchIndex>=0 && ((int)this.edges[matchIndex])==key) {
      leftIndex = matchIndex;
      matchIndex--;
    }
    
    ArrayList<Integer> retValues = new ArrayList<Integer>();
    while (leftIndex<getEdgeCount() && ((int)this.edges[leftIndex])==key) { 
      retValues.add((int)this.edges[leftIndex]);
      leftIndex++;
    }
    return retValues;
  }
    

  @Override
  public int getEdgeCount() {
    return edges.length;
  }


  @Override
  public EdgeRuler getEdgeRuler() {
    // just return the comparator defined there
    return Node.compositionEdgeRuler;
  }
}
