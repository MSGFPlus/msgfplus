package suffixtree.trees;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import msutil.AminoAcid;

import sequences.FastaSequence;
import sequences.ProteinFastaSequence;
import suffixtree.Constants;
import suffixtree.matches.ExactMatchObject;
import suffixtree.trees.deprecated.IntegerGappedSuffixTree;


/**
 * This data structure splits the sequence by initial letters, building a hash,
 * allowing for fast look up of the items and increasing the amount of amino
 * acids that can be stored by the data structure as a whole.
 * @author jung
 *
 */
public class HashedIntegerGappedSuffixTree {

  //private int hashSize;
  private FastaSequence sequence;
  HashMap<Integer,IntegerGappedSuffixTree> entries;
  HashMap<Integer,HashSet<Integer>> table;
  
  public void generate(HashMap<Integer,HashSet<Integer>> lookupTable) {
    Collection<Character> chars = this.sequence.getAlphabet();
    //System.err.println("Chars size: " + chars.size());
    //for (char c : chars) System.err.print(c + "");
    //System.err.println();
    HashMap<Integer,HashSet<String>> combos = new HashMap<Integer,HashSet<String>>();
    generateMassBins(chars, new ArrayList<Character>(), 0, combos);
    
    for (int mass : combos.keySet()) {
      HashSet<Integer> keys = new HashSet<Integer>();
      for(String s : combos.get(mass)) {
        
        
        if (s.length()>1) {
          ArrayList<Byte> sequence = new ArrayList<Byte>();
          for (int i=0; i<s.length(); i++) sequence.add(this.sequence.toByte(s.charAt(i)));
          chooseTwo(sequence, keys);
        }
        else {
          multiply(this.sequence.toByte(s.charAt(0)), this.sequence.getAlphabetAsBytes(), keys);
        }
      }
      //System.err.println(mass + " look up items: " + keys.size());
      lookupTable.put(mass, keys);
    }
  }
  
  
  private static void multiply(byte first, Set<Byte> alpha, HashSet<Integer> combos) {
    for (byte second : alpha) {
      combos.add((first<<8)|second);
    }
  }
  
  private static void chooseTwo(ArrayList<Byte> target, HashSet<Integer> combos) {
    for (int i=0; i<target.size(); i++) {
      for (int j=0; j<target.size(); j++) {
        if (i!=j) {
          combos.add((target.get(i)<<8)|target.get(j));
        }
      }
    }
  }
  
  
  private static void generateMassBins(Collection<Character> alpha, ArrayList<Character> prefix, int cumMass, HashMap<Integer,HashSet<String>> combos) {
  
    
    // add the current configuration to the combos
    StringBuffer sb = new StringBuffer();
    for (char letter : prefix) sb.append(letter);
    
    for (char letter : alpha) {
      ArrayList<Character> nextPrefix = new ArrayList<Character>(prefix);
      nextPrefix.add(letter);
      int nextCumMass = cumMass + AminoAcid.getStandardAminoAcid(letter).getNominalMass();
      if (nextCumMass<=Constants.MAX_GAP_MASS) 
        // recurse
        generateMassBins(alpha, nextPrefix, nextCumMass, combos);
    }
    
    if (cumMass>0) {
      if (!combos.containsKey(cumMass)) combos.put(cumMass,new HashSet<String>());
      combos.get(cumMass).add(sb.toString());
    }
  }
  
  
  
  /**
   * Build a HashedGappedSuffixTree using a hash size of 2 characters.
   * @param sequence
   */
  public HashedIntegerGappedSuffixTree(ProteinFastaSequence sequence) {
    this(sequence, 2);
  }

  /**
   * Constructor taking the sequence and the size of the hash.
   * @param sequence the sequence object to build this tree on.
   * @param hashSize the number of characters to do the hash on.
   */
  public HashedIntegerGappedSuffixTree(ProteinFastaSequence sequence, int hashSize) {
    
    this.sequence = sequence;
    
    this.table = new  HashMap<Integer,HashSet<Integer>>();
    generate(this.table);
    
    this.entries = new HashMap<Integer,IntegerGappedSuffixTree>();
    HashMap<Integer,Integer> counts = new HashMap<Integer,Integer>();
    
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
      
      int cumMass = sequence.getIntegerMass(start);
      for (int i=start+1; i<sequence.getSize(); i++) {
        if (this.sequence.isTerminator(i) || cumMass>Constants.MAX_QUERY_MASS) {
          break;
        }
        cumMass += sequence.getIntegerMass(i);
      }
      
      // too small to be inserted
      if (cumMass < Constants.MIN_QUERY_MASS) continue;
      
      if (!entries.containsKey(key)) {
        entries.put(key, new IntegerGappedSuffixTree(sequence, true));
        counts.put(key, 0);
      }
      
      counts.put(key, counts.get(key)+1);
      entries.get(key).insert(start);
    }
    
    System.out.println("Done inserting all suffixes to the trees");
    System.out.println(entries.size() + " entries in the HashedGappedSuffixTree");
    
    int count = 0;
    for (int key : entries.keySet()) {
      byte[] byteArray = new byte[hashSize];
      int target = key;
      for (int i=hashSize-1; i>=0; i--) {
        byteArray[i] = (byte)(target & 0xFF);
        target = target>>>8;
      }
      System.out.println(++count + ": Processing sequences starting with: " + this.sequence.toString(byteArray) + " " + counts.get(key));
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
  public void search(ArrayList<Integer> query, HashSet<ExactMatchObject> results) {
    
    int firstMass = query.get(0);
    if (this.table.containsKey(firstMass)) {
      for(int key : this.table.get(firstMass)) {
        if (this.entries.containsKey(key)) {
          this.entries.get(key).search(query, results);
        }
      }
    }
    
    return;
  }
  
}
