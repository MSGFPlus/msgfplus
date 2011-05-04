package suffixtree.trees;

import java.util.HashSet;

import msutil.AminoAcid;

import sequences.ProteinFastaSequence;
import suffixtree.Constants;
import suffixtree.edges.Edge;
import suffixtree.edges.ByteEdge;
import suffixtree.nodes.InternalNode;
import suffixtree.nodes.Node;


/**
 * This is the data structure for storing the database as integer masses.
 * @author jung
 *
 */
public class IntegerSuffixTree {

  
  
  /**
   * Helper class that defines an edge in the graph taking two integers for the
   * start and end indices of the sequence. The reason this is an inner class 
   * is because this edges needs access to the backing FastaSequence.
   * @author jung
   *
   */
  public class CompressedEdge extends Edge {
    
    private int start;
    private int end;
    
    /**
     * Constructor.
     * @param sink the destination node
     * @param start the start index in the FastaSequence.
     * @param end the end position in the FastaSequence.
     */
    public CompressedEdge(Node sink, int start, int end) {
      //System.out.println("Created " + sequence.toS)
      this.start = start;
      this.end = end;
      setSink(sink);
    }

    @Override
    public int getLabelAt(int offset) {
      return sequence.getIntegerMass(start+offset);
    }

    @Override
    public int size() {
      return end-start;
    }

    @Override
    public Edge split(int offset) {
      if (this.start+offset >= this.end) return null;
      Edge otherEdge = new CompressedEdge(getSink(), start+offset, end);
      this.end = this.start + offset;
      
      setSink(new InternalNode(otherEdge));
      return otherEdge;
    }
    
    @Override
    public String toString() {
      StringBuffer result = new StringBuffer();
      result.append("IntegerCompressedEdge: " + sequence.getSubsequence(start, Math.min(start+6,end)));
      result.append("[" + start + "," + end + "). Positions:");
      for (int pos : getSink().getPositions()) {
        result.append(" " + pos);  
      }
      result.append(" #" + hashCode());
      return result.toString();
    }

    @Override
    public int getEnd() {
      return end;
    }
    
    @Override
    public int getStart() {
      return start;
    }

    @Override
    public int length() {
      return size();
    }
    
    @Override
    public int compareTo(Edge o) {
      // Lexicographical ordering
      int minLength = Math.min(this.size(), o.size());
      for (int i=0; i<minLength; i++) {
        if(this.getLabelAt(i) < o.getLabelAt(i)) return -1;
        if(this.getLabelAt(i) > o.getLabelAt(i)) return 1;
      }
      
      if (this.size() < o.size()) return -1;
      if (this.size() > o.size()) return 1;
      return 0;
    }

    @Override
    public int mass() {
      String slice = sequence.getSubsequence(start, end);
      //System.out.println(slice);
      int cumMass = 0;
      for (int i=0; i<slice.length(); i++) {
        //System.out.println(slice.charAt(i));
        cumMass += AminoAcid.getStandardAminoAcid(slice.charAt(i)).getNominalMass();
      }
      return cumMass;
    }
    
    @Override
    public CompressedEdge getClone() {
      return new CompressedEdge(getSink(), this.start, this.end);
    }
  }
  
  
  
  // Members
  private ProteinFastaSequence sequence;      // the backing sequence
  private InternalNode root;           // the top of the tree
  
  /**
   * Constructor taking a ProteinFastaSequence object.
   * @param sequence the sequence to build the graph on.
   */
  public IntegerSuffixTree(ProteinFastaSequence sequence) {
    this.sequence = sequence;  
    this.root = null;
    for (int start=0; start<sequence.getSize(); start++) insert(start);
  }
  
  /**
   * Special constructor creating an empty tree and not inserting anything yet.
   * This allows to selective insertion later on by calling the insert method.
   * @param sequence the FastaSequence object.
   * @param noInsert dummy parameter to distinguish this constructor from the 
   *                 other constructor.
   */
  public IntegerSuffixTree(ProteinFastaSequence sequence, boolean noInsert) {
    this.sequence = sequence;
    this.root = null;
  }
  
  @Override
  public String toString() {
    StringBuffer sb = new StringBuffer();
    for (int i=0; i<root.getDegree(); i++) {
      sb.append("Edge: " + root.getEdgeAt(i).getLabel() + "\n");
    }
    return sb.toString(); 
  }
  
  /**
   * Insert a given position of the stored FastaSequence into the tree. 
   * @param start
   */
  protected void insert(int start) {
    // don't insert invalid position
    if (start>=this.sequence.getSize() || start<0) return;
    
    if (!this.sequence.hasMass(start)) return;
    
    int end = start;
    int cumMass = this.sequence.getIntegerMass(start);
    for (int i=start+1; i < sequence.getSize(); i++) {
      if (!this.sequence.hasMass(i) || cumMass>Constants.MAX_QUERY_MASS) {
        end = i;
        break;
      }
      cumMass += this.sequence.getIntegerMass(i);
    }
    
    // check that the sequence falls within the range of allowable masses
    if (cumMass<Constants.MIN_QUERY_MASS) return;
    
    if (this.root==null) {
      // initialize root
      this.root = new InternalNode(new CompressedEdge(new InternalNode(start), start, end));
    }
    else {
      Edge currentEdge = new CompressedEdge(new InternalNode(start), start, end);
      root.insert(currentEdge);
    }
  }
  
  /**
   * Getter method to retrieve the root.
   * @return the root of this tree.
   */
  public InternalNode getRoot() {
    return root;
  }
  
  /**
   * Getter method to retrieve the sequence backing this tree. 
   * @return the ProteinFastaSequencing backing this tree.
   */
  public ProteinFastaSequence getSequence() {
    return this.sequence;
  }
  
  
  public HashSet<Integer> search(byte[] query) {
    Node current = root;
    
    System.err.println("Broken Search function");
    HashSet<Integer> results = new HashSet<Integer>();
    
    for (byte item : query) {
      Edge queryEdge = new ByteEdge(item);
      int match = current.search(queryEdge);
      if (match<0) {
        return null;
      }
      Edge path = current.getEdgeAt(match);
      if (path.size() < queryEdge.size()) {
        queryEdge = queryEdge.split(path.size());
        current = path.getSink();
      }
      else {
        path.getSink().getAllPositions(results);
        break;
      }
    }
    return results;
  }
  
  
  public CompressedEdge createCompressedEdge(InternalNode sink, int start, int end) {
    return new CompressedEdge(sink, start, end);
  }
  
}
