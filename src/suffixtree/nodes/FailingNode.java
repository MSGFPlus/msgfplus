package suffixtree.nodes;

public class FailingNode extends InternalNode {
  private InternalNode fail;
  
  public FailingNode(int position) {
    super(position);
  }
  
  public FailingNode() {
    super();
  }
  
  public FailingNode(InternalNode n) {
    super(n.getEdges(), n.getDegree(), n.getPositions());
    if (n instanceof FailingNode) {
      this.fail = ((FailingNode)n).fail;
    }
  }
  
  public Node getFail() {
    return fail;
  }
  
  public void setFail(InternalNode n) {
    this.fail = n;
  }

}
