package suffixtree.edges;


/**
 * This Object is to create edges for querying a graph.
 * @author jung
 *
 */
public class QueryEdge extends Edge {

  int label;
  
  /**
   * Constructor taking a byte array.
   * @param sequence the byte array as the label of this edge.
   */
  public QueryEdge(int label) {
    this.label = label;
  }
  
  @Override
  public int getLabelAt(int offset) {
    return label;
  }

  @Override
  public int size() {
    return 1;  
  }
  
  @Override
  public int length() {
    return 1;
  }

  @Override
  public Edge split(int offset) {
    System.err.println("Unsupported operation: Split query edge");
    return null;
  }

  @Override
  public int getEnd() {
    System.err.println("Unsupported operation");
    System.exit(-1);
    return 0;
  }

  @Override
  public int getStart() {
    System.err.println("Unsupported operation");
    System.exit(-1);
    return 0;
  }

  @Override
  public int mass() {
    return label;
  }
}
