package suffixtree.edges;

import java.util.ArrayList;
import java.util.List;

import suffixtree.nodes.ComplexInternalNode;
import suffixtree.nodes.Node;

/**
 * This class represents a query edge backed up by an array of integers 
 * representing the labels of the edge. Note that this edge is still a single
 * object, but might have multiple labels. This edge can be split into smaller
 * portions when building the tree by calling the split method.
 * @author jung
 *
 */
public class DirectedMassEdge extends MassEdge {

  public DirectedMassEdge(ArrayList<Integer> labels, Node sink) {
    super(labels, sink);
  }
    
  @Override 
  public DirectedMassEdge duplicate() {
    return new DirectedMassEdge(new ArrayList<Integer>(getLabels()), getSink());
  }
 
  @Override
  public Edge split(int offset) {
    if (offset >= this.size()) return null;
    
    ComplexInternalNode thisSink = (ComplexInternalNode)getSink();
    List<Integer> tempList = this.getLabels().subList(offset, this.getLabels().size());
    DirectedMassEdge second = new DirectedMassEdge(new ArrayList<Integer>(tempList), thisSink);
    tempList.clear(); // remove the items from the first list

    // the new middle node
    ComplexInternalNode middle = new ComplexInternalNode(second);
    middle.setParentNode(thisSink.getParentNode());
    setSink(middle);
    thisSink.setParentNode(middle);

    return second;
  }

  
  /**
   * Retrieves the 0-based index of the edge that meets this mass.
   * @param mass the mass to meet
   * @return the index of the edge that meets the mass from the start of the edge array.
   *         for example if the edge array is [57, 71, 57], an input parameter 57 will
   *         return 0.
   */
  public int getIndexFromStart(int mass) {
    int cumMass = 0;
    for (int i=0; i<getLabels().size(); i++) {
      cumMass += getLabelAt(i);
      if (cumMass>=mass) return i;
    }
    return getLabels().size();
  }
  
  
  /**
   * Retrieves the 0-based index of the edge that meets this mass, accumulating 
   * from the end (last edge).
   * @param mass the mass to meet
   * @return the index of the edge that meets the mass from the end of the edge array.
   *         for example if the edge array is [71, 71, 57], an input parameter 57 will
   *         return 2.
   */
  public int getIndexFromEnd(int mass) {
    int cumMass = 0;
    for (int i=getLabels().size()-1; i>=0; i--) {
      cumMass += getLabelAt(i);
      if (cumMass>=mass) return i;
    }
    return 0;
  }
  
  
  /**
   * Retrieve the cumulative mass at the given index, accumulating from the beginning.
   * @param index the index to stop accumulating (inclusive).
   * @return the mass accumulated.
   */
  public int getMassAtIndexForward(int index) {
    int cumMass = 0;
    for (int i=0; i<=index; i++) {
      cumMass += getLabelAt(i);
    }
    return cumMass;
  }
  
  
  /**
   * Retrieve the cumulative mass at the given index, accumulating from the end.
   * @param index the index to stop accumulating (inclusive).
   * @return the mass accumulated
   */
  public int getMassAtIndexBackward(int index) {
    int cumMass = 0;
    for (int i=getLabels().size()-1; i>=index; i--) {
      cumMass += getLabelAt(i);
    }
    return cumMass;
  }
}
