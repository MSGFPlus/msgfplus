package suffixtree.edges;

import java.util.ArrayList;
import java.util.List;

import suffixtree.nodes.FailingNode;
import suffixtree.nodes.InternalNode;
import suffixtree.nodes.Node;

public class MassEdge extends Edge {

  private ArrayList<Integer> labels;
  
  public MassEdge(ArrayList<Integer> labels, Node sink) {
    this.labels = new ArrayList<Integer>(labels);
    setSink(sink);
  }
  
  public MassEdge(ArrayList<Integer> labels) {
    this.labels = new ArrayList<Integer>(labels);
  }
  
  public MassEdge duplicate() {
    return new MassEdge(new ArrayList<Integer>(this.labels), getSink());  
  }
  
  public MassEdge(int[] labels, InternalNode sink) {
    this.labels = new ArrayList<Integer>();
    for (int label : labels) this.labels.add(label);
    setSink(sink);
  }
  
  /**
   * Getter method that gets the list of masses of this edge.
   * @return the masses of this object as an array
   */
  public ArrayList<Integer> getLabels() {
    return labels;
  }
  
  /**
   * Calculate the total mass of this edge
   * @return the total mass of this edge
   */
  public int getTotalMass() {
    int result = 0;
    for (int mass : this.labels) result += mass;
    return result;
  }
  
  /**
   * This converts the Tree to a Trie.
   */
  public void divideAll() {
    if (this.size()==1) return;
    
    // create individual edges for each mass
    FailingNode currentNode = null;
    MassEdge edge = null;
    MassEdge subEdge = null;
    for (int mass : this.labels) {
      int[] massArray = {mass};
      FailingNode sinkNode = new FailingNode();
      subEdge = new MassEdge(massArray, sinkNode);
      if (currentNode!=null) currentNode.insert(subEdge, 0);
      else                   edge = subEdge;
      currentNode = sinkNode;
    }
    this.labels = edge.labels;
    subEdge.setSink(this.getSink());
    this.setSink(edge.getSink());
    
  }
  
  @Override
  public int getEnd() {
    System.err.println("IntegerEdge does not support getting the end position");
    return 0;
  }

  @Override
  public int getLabelAt(int offset) {
    return labels.get(offset);
  }

  @Override
  public int getStart() {
    System.err.println("IntegerEdge does not support getting the start position");
    return 0;
  }

  @Override
  public int length() {
    return labels.size();
  }

  @Override
  public int size() {
    return labels.size();
  }

  @Override
  public Edge split(int offset) {
    if (offset >= this.size()) return null;
    List<Integer> tempList = this.labels.subList(offset, this.labels.size());
    MassEdge second = new MassEdge(new ArrayList<Integer>(tempList), getSink());
    tempList.clear(); // remove the items from the first list
    setSink(new InternalNode(second));
    return second;
  }

  @Override
  public int mass() {
    System.err.println("MassEdge does not support mass(). Ironic huh?");
    return 0;
  }
  
  @Override
  public String toString() {
    StringBuffer sb = new StringBuffer();
    sb.append("MassEdge: {");
    for (int mass : this.labels) {
      sb.append(mass + ",");
    }
    sb.deleteCharAt(sb.length()-1);
    sb.append("}");
    return sb.toString();
    
  }

}
