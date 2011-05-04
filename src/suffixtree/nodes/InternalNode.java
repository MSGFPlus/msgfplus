package suffixtree.nodes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Stack;

import suffixtree.edges.Edge;


/**
 * The generalized data structure representing a node inside the KeywordTree.  
 * @author jung
 *
 */
public class InternalNode implements Node {

  // all the outgoing edges from here
  private Edge[] edges;
  private int degree;
  
  // all the starting positions stored by this node
  private int[] positions;
  
  
  public InternalNode duplicate() {
    InternalNode result = new InternalNode();
    result.edges = new Edge[this.edges.length];
    for (int i=0; i<result.edges.length; i++) result.edges[i] = this.edges[i];
    result.positions = Arrays.copyOf(this.positions, this.positions.length);
    result.degree = this.degree;
    return result;
  }
  
  
  public InternalNode(ArrayList<Integer> positions) {
    this.edges = new Edge[1];
    this.degree = 0;
    this.positions = new int[positions.size()];
    for (int i=0; i<this.positions.length; i++) this.positions[i] = positions.get(i);
  }
  
  public InternalNode() {
    this.edges = new Edge[1];
    this.degree = 0;
    this.positions = null;
  }
  
  public InternalNode(Edge[] edges, int degree, int[] positions) {
    this.edges = edges;
    this.degree = degree;
    this.positions = positions;
  }
  
  public InternalNode(Edge e) {
    this.edges = new Edge[1];
    this.edges[0] = e;
    this.degree = 1;
    this.positions = null;
  }
  
  /**
  public InternalNode(Edge[] edges, int degree) {
    this.edges = edges;
    this.degree = degree;
    this.positions = new int[0];
  } **/
  
  public InternalNode(int position) {
    this.edges = new Edge[1];
    this.degree = 0;
    this.positions = new int[1];
    this.positions[0] = position;
  }
  
  /**
   * Standard getter method.
   * @return the array of edges stores by this node. Note that not all positions
   *         of the array are used. To find out use the getDegree method.
   */
  public Edge[] getEdges() {
    return this.edges;
  }
  
  @Override
  public void insert(Edge edge) {
    // find the insertion point
    int pos = search(edge);
    
    if (pos<0) {
      // not found just insert this new edge into the edge array
      insert(edge, -pos-1);
    }
    else {
      // found and the matching edge is at pos
      int lcp = edges[pos].getLongestCommonPrefixSize(edge);
      
      // this split operation will rearrange the tree properly
      edges[pos].split(lcp);
      
      // now split the inserting edge if applicable
      Edge otherEdge = edge.split(lcp);
      if (otherEdge==null) {
        // simply add end positions to the sink of the split existing edge
        for (int x : edge.getSink().getPositions()) edges[pos].getSink().addPosition(x);
      }
      else {
        // insert the split edge into the new internal node
        edges[pos].getSink().insert(otherEdge);
      }
    }
  }
  
  @Override
  public void insert(Edge edge, int index) {
    if (this.degree>=this.edges.length) {
      // need to double the capacity of this set of edges
      Edge[] tempEdges = new Edge[Math.max(this.edges.length*2, 1)];
      System.arraycopy(this.edges, 0, tempEdges, 0, index);
      System.arraycopy(this.edges, index, tempEdges, index+1, this.degree-index);
      tempEdges[index] = edge;
      this.edges = tempEdges;
    }
    else {
      System.arraycopy(this.edges, index, this.edges, index+1, this.degree-index);  
      this.edges[index] = edge;
    }
    this.degree++;
  }

  @Override
  public int search(Edge edge) {
    // we know all edges must have different starting positions, so we only need
    // to look at the start position
    int lower = 0, upper = this.degree;
    
    while (true) {
      // base case, not found
      if (lower==upper) return -lower-1;
      
      // base case, possibly match
      if (upper-lower==1) {
        int targetChar = this.edges[lower].getLabel(); 
        if (targetChar < edge.getLabel()) {
          // not found, insert at the end position
          return -upper-1;
        }
        else if (targetChar > edge.getLabel()) {
          // not found, insert at the start position
          return -lower-1;
        }
        else {
          // found
          return lower;
        }
      }
      
      // recursive search
      int mid = (lower+upper) / 2;
      int midByte = this.edges[mid].getLabel();
      if (edge.getLabel() < midByte) {
        // recurse on the top half
        upper = mid;
      }
      else if (edge.getLabel() > midByte) {
        // recurse of the bottom half
        lower = mid+1;
      }
      else {
        // found
        return mid;
      }
    }
  }
  
  @Override
  public int search(int key, int lower, int upper) {
    while (true) {
      // base case, not found
      if (lower==upper) return -lower-1;
      
      // base case, possibly match
      if (upper-lower==1) {
        int targetChar = this.edges[lower].getLabel(); 
        if (targetChar < key) {
          // not found, insert at the end position
          return -upper-1;
        }
        else if (targetChar > key) {
          // not found, insert at the start position
          return -lower-1;
        }
        else {
          // found
          return lower;
        }
      }
      
      // recursive search
      int mid = (lower+upper) / 2;
      int midByte = this.edges[mid].getLabel();
      if (key < midByte) {
        // narrow the upper limit
        upper = mid;
      }
      else if (key > midByte) {
        // recurse of the bottom half
        lower = mid+1;
      }
      else {
        // found
        return mid;
      }
    }
  }
  
  @Override
  public int search(int key, int lower) {
	  return search(key, lower, this.degree);
  }
  
  @Override
  public int search(int key) {
	  return search(key, 0, this.degree);
  }

  @Override
  public void addPosition(int position) {
    
    if (this.positions==null) {
      this.positions = new int[1];
      this.positions[0] = position;
      return;
    }
    
    for (int pos : this.positions) {
      if (pos==position) return;
    }
    
    this.positions = Arrays.copyOf(this.positions, this.positions.length+1);
    this.positions[positions.length-1] = position;
  }
  
  @Override
  public void addPositions(ArrayList<Integer> positions) {
    for (int pos : positions) addPosition(pos);
  }
  
  @Override
  public int[] getPositions() {
    if (positions==null) return new int[0];
    return positions;
  }

  @Override
  public int getPositionsCount() {
    if (positions==null) return 0;
    return positions.length;
  }
  
  @Override
  public Node setEdges(Edge[] edges) { 
  	this.edges = edges;
  	this.degree = edges.length;
  	return this;
  }
  
  @Override
  public Node setEdges(Edge[] edges, int degree) { 
    this.edges = edges;
    this.degree = degree;
    return this;
  }

  @Override
  public Edge getEdgeAt(int index) {
    return edges[index];
  }

  @Override
  public void setEdgeAt(int index, Edge e) {
    edges[index] = e;  
  }
  
  @Override
  public void getAllPositions(HashSet<Integer> results) {
    HashSet<Node> seen = new HashSet<Node>();
    Stack<Node> nodes = new Stack<Node>();
    nodes.push(this);
    
    while (!nodes.isEmpty()) {
      Node nextNode = nodes.pop();
      if (!seen.contains(nextNode)) {
        seen.add(nextNode);
        for(int position : nextNode.getPositions()) {
          results.add(position);
        }
        // add the children
        for (int i=0; i<nextNode.getDegree(); i++) {
          nodes.add(nextNode.getEdgeAt(i).getSink());
        }
      }
    }
    return;
  }

  
  @Override
  public int getDegree() {
    return this.degree;
  }
  
  @Override
  public String toString() {
    String result = "InternalNode[" + degree + "] with start positions: ";
    for (int i : this.positions) {
      result += i + " ";
    }
    return result;
  }


  @Override
  public void setPositions(int[] positions) {
    this.positions = positions;
  }
  

  @Override
  public Edge getMaximumEdge() {
    if (this.getDegree()>0) return this.getEdgeAt(this.getDegree()-1);
    return null;
  }


  @Override
  public Edge getMinimumEdge() {
    if (this.getDegree()>0) return this.getEdgeAt(0);
    return null;
  }

}
