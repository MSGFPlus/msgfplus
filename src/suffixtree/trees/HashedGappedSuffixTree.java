package suffixtree.trees;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import sequences.FastaSequence;
import suffixtree.Constants;
import suffixtree.edges.ByteEdge;



/**
 * This is a easy memory saver that scales the GappedSuffixTree to databases of
 * sizes up to 5 million with machines with 8GB of RAM with minimum impact on
 * the query time (slow down of less than 10X).
 * @author jung
 *
 */
public class HashedGappedSuffixTree {

  //private int hashSize;
  private FastaSequence sequence;
  private HashMap<Integer,GappedSuffixTree> entries;
  
  /**
   * Build a HashedGappedSuffixTree using a hash size of 2 characters.
   * @param sequence
   */
  public HashedGappedSuffixTree(FastaSequence sequence) {
    this(sequence, 2);
  }

  
  
  /**
   * Constructor taking the sequence and the size of the hash.
   * @param sequence the sequence object to build this tree on.
   * @param hashSize the number of characters to do the hash on.
   */
  public HashedGappedSuffixTree(FastaSequence sequence, int hashSize) {
    
    //this.hashSize = hashSize;
    this.sequence = sequence;
    this.entries = new HashMap<Integer,GappedSuffixTree>();
    HashMap<Integer,ArrayList<Integer>> positions = new HashMap<Integer,ArrayList<Integer>>();
    
    if (hashSize > 4) {
      System.err.println("hashSize of " + hashSize + " is currently not supported.");
      System.exit(-1);
    }
    
    for (int start=0; start<sequence.getSize()-hashSize; start++) {
      int key = 0;
      boolean completeFlag = true;
      for (int i=0; i<hashSize; i++) {
        byte b = sequence.getByteAt(start+i);
        if (sequence.isTerminator(start+i)) {
          completeFlag = false;
          break;
        }
        key = (key<<8) | b;
      }
      
      if (!completeFlag)     continue;
      
      int end = start;
      for (int i=start+1; i < sequence.getSize(); i++) {
        if (this.sequence.isTerminator(i)) {
          end = i;
          break;
        }
      }
      
      // too small to be inserted
      if (end-start < Constants.MIN_QUERY_CHAR) continue;
      
      if (!entries.containsKey(key)) {
        entries.put(key, new GappedSuffixTree(sequence, true));
        positions.put(key, new ArrayList<Integer>());
      }
      entries.get(key).insert(start);
      positions.get(key).add(start);
    }
    
    System.out.println("Done inserting all suffixes to the trees");
    System.out.println(entries.size() + " entries in the HashedGappedSuffixTree");
    
    for (int key : entries.keySet()) {
      byte[] byteArray = new byte[hashSize];
      int target = key;
      for (int i=hashSize-1; i>=0; i--) {
        byteArray[i] = (byte)(target & 0xFF);
        target = target>>>8;
      }
      System.out.print("Processing sequences starting with: " + this.sequence.toString(byteArray) + " " + key + " ");
      System.out.println(positions.get(key).size() + " items");
      entries.get(key).makeGapLinks();
    }
  }
  
  
  
  /**
   * Getter method for the sequence in which this tree is based on.
   * @return the FastaSequence object.
   */
  public FastaSequence getSequence() {
    return sequence;
  }
  

  
  /**
   * Search for a query
   * @param query the query object.
   * @param results store the results here.
   */
  public void search(ArrayList<ByteEdge> query, HashSet<Integer> results) {
    // lazy implementation, just query everything    
    for(GappedSuffixTree t : entries.values()) {
      t.search(query, results);
    }
    return;
  }
  
}
