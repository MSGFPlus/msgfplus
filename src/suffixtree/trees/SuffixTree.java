package suffixtree.trees;

import java.util.HashSet;

import sequences.FastaSequence;
import suffixtree.Constants;
import suffixtree.edges.Edge;
import suffixtree.edges.QueryEdge;
import suffixtree.nodes.InternalNode;
import suffixtree.nodes.Node;


/**
 * General suffix tree that takes as input a FastaSequence.
 * @author jung
 *
 */
public class SuffixTree {

  private FastaSequence sequence;
  private InternalNode root;
  
  /**
   * Helper class that defines an edge in the graph taking two integers for the
   * start and end indices of the sequence.
   * @author jung
   *
   */
  public class CompressedEdge extends Edge {
    
    private int start;
    private int end;
    
    public CompressedEdge(Node sink, int start, int end) {
      this.start = start;
      this.end = end;
      setSink(sink);
    }

    @Override
    public int getLabelAt(int offset) {
      return sequence.getByteAt(start+offset);
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
      return "CompressedEdge: " + sequence.getSubsequence(start, Math.min(start+6,end)) + "[" + start + "," + end + ")" + " #" + hashCode();  
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
    public int mass() {
      System.err.println("Unsupported method: SuffixTree.CompressedEdge.mass()");
      return 0;
    }
  }
  

  /**
   * Constructor taking a FastaSequence object.
   * @param sequence the sequence to build the graph on.
   */
  public SuffixTree(FastaSequence sequence) {
    this.sequence = sequence;  
    this.root = null;
    for (int start=0; start<sequence.getSize(); start++) {
      insert(start);
    }
  }
  
  
  /**
   * Special constructor creating an empty tree and not inserting anything yet.
   * This allows to selective insertion later on by calling the insert method.
   * @param sequence the FastaSequence object.
   * @param noInsert dummy parameter to distinguish this constructor from the 
   *                 other constructor.
   */
  public SuffixTree(FastaSequence sequence, boolean noInsert) {
    this.sequence = sequence;
    this.root = null;
  }
  
  
  public void insert(int start) {
    if (start>=this.sequence.getSize()) return;
    
    if (this.sequence.isTerminator(start)) return;
    
    int end = start;
    for (int i=start+1; i<sequence.getSize(); i++) {
      if (this.sequence.isTerminator(i)) {
        end = i;
        break;
      }
    }
    end = Math.min(start+Constants.MAX_QUERY_CHAR, end);
    
    // too small to be inserted
    if (end-start < Constants.MIN_QUERY_CHAR) return;
    
    if (this.root==null) {
      this.root = new InternalNode(new CompressedEdge(new InternalNode(start), start, end));
    }
    else {
      Edge currentEdge = new CompressedEdge(new InternalNode(start), start, end);
      root.insert(currentEdge);
    }
  }
  
  public InternalNode getRoot() {
    return root;
  }
  
  public FastaSequence getSequence() {
    return this.sequence;
  }
  
  
  public void search(byte[] query, HashSet<Integer> results) {
    Node current = root;
    
    int i = 0;
    while (true) {
      int matchIndex = current.search(new QueryEdge(query[i])); 
      if (matchIndex < 0) return;
      
      i++;
      if (i==query.length) {
        current.getAllPositions(results);
        return;
      }
      current = current.getEdgeAt(matchIndex).getSink();
    }
  }
  
  
  /*
  public void searchExact(MassEdge query, int index, InternalNode target, HashSet<Integer> matches) {
    if (target==null) target = this.root;
    
    // base case
    if (index >= query.size()) {
      //System.out.println("Done at index " + index);
      target.getAllPositions(matches);
      return;
    }
    
    //System.out.println(index + " Querying " + query.getLabelAt(index));
    
    int matchIndex = target.search(new QueryEdge(query.getLabelAt(index)));
    if (matchIndex < 0) {
      // no match, try next mass
      System.out.println("Not found " + query.getLabelAt(index));
      System.out.println("Outgoing edges of mis matching node:");
      for (int i=0; i<target.getDegree(); i++) {
        System.out.println(target.getEdgeAt(i));
      }
      System.out.println();
      return;
    }
    
    // advance the index greedily until you reach the end of the matching edge
    Edge match = target.getEdgeAt(matchIndex);
    int advance = 1;
    for (; advance<match.size(); advance++) {
      
      // done matching
      if (advance+index>=query.length()) return;
      
      // no match
      if (match.getLabelAt(advance)!=query.getLabelAt(advance+index)) return; 
    }
    
    // recurse
    searchExact(query, index+advance, target.getEdgeAt(matchIndex).getSink(), matches);  
  }
  */

  
  
  public CompressedEdge createCompressedEdge(InternalNode sink, int start, int end) {
    return new CompressedEdge(sink, start, end);
  }
  
  
  public static void main(String[] args) {
    String fastaFile;
    String userHome = System.getProperty("user.home");
    
    fastaFile = userHome+"/Data/Databases/test.fasta";
    fastaFile = userHome+"/Data/Databases/small.fasta";
    //fastaFile = userHome+"/Data/Databases/medium.fasta";
    //fastaFile = userHome+"/Data/Databases/million.fasta";
    //fastaFile = userHome+"/Data/Databases/large.fasta";
    //fastaFile = userHome+"/Data/Databases/uniprot_sprot.fasta";
    //fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    //fastaFile = userHome+"/Data/Databases/ShewDB/SOne_proteins_withContams.fasta";
    
    //String graphFile = fastaFile.replaceAll(".fasta$", "")+".st";
    
    long time = System.currentTimeMillis();
    FastaSequence sequence = new FastaSequence(fastaFile, sequences.Constants.AMINO_ACIDS_18);
    System.out.println("-- Loading fasta file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    System.out.println("--- Number of characters in fasta file: " + sequence.getSize());
    
    time = System.currentTimeMillis();
    SuffixTree st = new SuffixTree(sequence);
    System.out.println("-- Loading SuffixTree file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    //System.out.println(st);
    
    // Query all possible sequences
    for (int start=1; start < sequence.getSize(); start++) {
      
      if (sequence.isTerminator(start)) continue;
      
      // for each start position, find the end
      int end = start;
      for (int i=start+1; i < sequence.getSize(); i++) {
         if (sequence.isTerminator(i)) {
          end = i;
          break;
        }
      }
      
      if (end > start+10) {
        
        HashSet<Integer> matches = new HashSet<Integer>();
        st.search(sequence.getBytes(start, end), matches);
        if (matches.size()!=0) {
          //System.err.println(start + "\t" + matches.iterator().next());
        }
        else {
          System.err.println("Not found " + sequence.toString(sequence.getBytes(start, end)));
        }
      }
      
      break;
    }
  }
  
}
