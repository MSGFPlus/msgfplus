package suffixtree.nodes;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Stack;

import suffixtree.edges.Edge;

public class FinalInternalNode implements Node {
  
  private Edge[] edges;      // the edge array indicating whether it exists
  private int[] positions;   // number of positions
  private Edge minEdge;      // the minEdge label
  private Edge maxEdge;      // the maxEdge label 

  /**
   * Constructor that converts the given InternalNode into a final node that
   * does not allow modifications and has fast (constant) access to edges.
   * @param n the node in which this node is based on
   */
  public FinalInternalNode(Node n) {
    
    this.positions = n.getPositions();
    
    if (n.getDegree()>0) {
      this.edges = new Edge[n.getEdgeAt(n.getDegree()-1).getLabel()+1];
    }
    else {
      this.edges = new Edge[0]; // empty node
    }
    
    for (int i=0; i<n.getDegree(); i++) {
      this.edges[n.getEdgeAt(i).getLabel()] = n.getEdgeAt(i);
    }

    this.minEdge = n.getMinimumEdge();
    this.maxEdge = n.getMaximumEdge();
  }
  
  @Override
  public void addPosition(int position) {
    System.err.println("Cannot insert new position into this node");
    System.exit(-9);
  }

  @Override
  public void addPositions(ArrayList<Integer> positions) {
    System.err.println("Cannot insert new positions into this node");
    System.exit(-9);
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
    // this is pretty meaningless for this type of nodes
    return edges.length;
  }

  @Override
  public Edge getEdgeAt(int index) {
    return this.edges[index];
  }

  @Override
  public int[] getPositions() {
    if (this.positions==null) return new int[0];
    return this.positions;
  }

  @Override
  public int getPositionsCount() {
    if (this.positions==null) return 0;
    return this.positions.length;
  }
  
  @Override
  public void insert(Edge edge, int index) {
    System.err.println("Cannot insert new edge into this node");
    System.exit(-9);
  }

  @Override
  public void insert(Edge edge) {
    System.err.println("Cannot insert new edge into this node");
    System.exit(-9);
  }

  @Override
  public int search(Edge edge) {
    return search(edge.getLabel());
  }

  @Override
  public int search(int key) {
    if (key >= this.edges.length || this.edges[key]==null) return -key-1;
    return key;
  }

  @Override
  public int search(int key, int lower) {
    return search(key);
  }
  
  @Override
  public int search(int key, int lower, int upper) {
    return search(key);
  }
  
  @Override
  public void setEdgeAt(int index, Edge e) {
    this.edges[index] = e;
  }

  @Override
  public Node setEdges(Edge[] edges) {
    System.err.println("Cannot modify this node by replacing its edge array");
    System.exit(-9);
    return null;
  }

  @Override
  public Node setEdges(Edge[] edges, int degree) {
    System.err.println("Cannot modify this node by replacing its edge array");
    System.exit(-9);
    return null;
  }

  @Override
  public void setPositions(int[] positions) {
    System.err.println("Cannot modify this node by replacing its position array");
    System.exit(-9);
  }

  @Override
  public Edge getMaximumEdge() {
    return this.maxEdge;
  }

  @Override
  public Edge getMinimumEdge() {
    return this.minEdge;
  }

}
