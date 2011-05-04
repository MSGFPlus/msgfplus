package gappedsequences;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.TreeMap;

import sequences.ProteinFastaSequence;
import suffixtree.Constants;
import suffixtree.trees.KeywordTree;

/**
 * This is the k-mer represented by a starting position and a bit vector 
 * of how the splits.
 * @author jung
 *
 */
public class SequenceFactory {
  
  private class Sequence implements Comparable<Sequence>, Iterable<Integer> {
    
    /**
     * Allows the iteration of the this sequence, returning masses
     * @author jung
     *
     */
    private class MassIterator implements Iterator<Integer> {
      private int currentPosition;
      MassIterator() {
        currentPosition = 0;
      }
      
      @Override
      public boolean hasNext() {
        return currentPosition < length;
      }

      @Override
      public Integer next() {
        // count how many shifts before we encounter the left most 0 after currentPosition
        int cumMass = 0;
        int number = 0;
        do {
          cumMass += sequence.getIntegerMass(start+currentPosition);
          number = (breaks>>>(31-currentPosition++))&1;
          //System.out.println("number is " + number + " position is " + currentPosition);
        } while (number!=0);
        return cumMass;
      }

      @Override
      public void remove() {
        System.err.println("This iterator does not support this operation");
      }
    }
    
    private long start;
    private int length;
    private int breaks;
    private int maxMass;
    private int massCount;
    //private boolean unique;
    
    Sequence(long start, int length, int breaks) {
      this.start = start;
      this.length = length;
      this.breaks = breaks;
      //this.unique = false;
      
      Iterator<Integer> it = iterator();
      while (it.hasNext()) {
        this.massCount++;
        int mass = it.next();
        if (mass>this.maxMass) this.maxMass=mass;
      }
    }
    
    public ArrayList<Integer> toArray() {
      ArrayList<Integer> ret = new ArrayList<Integer>();
      Iterator<Integer> it = this.iterator();
      while (it.hasNext()) {
        ret.add(it.next());
      }
      return ret;
    }
    
    /**
     * Derive a new sequence from the current sequence, by merging masses
     * from the 'start'-index up to 'extend'-positions
     * @param start the index of the mass of this sequence to start merging
     * @param extend the extension count to merge masses
     * @return the new sequence
     */
    public Sequence derive(int start, int extend) {
      int breaks = this.breaks;
      int position = 0;
      
      // first move start positions
      int skip = 0;
      for(int index=31; index>=0; index--) {
        if (position>=start) break;
        if (((breaks>>>index)&1)==0) position++;
        skip++;
      }

      //System.out.println("Position " + skip);
      
      int count = 0;
      for(int index=31-skip; index>=0; index--) {
        if (((breaks>>>index)&1)==0) { 
          // set the current index to 1
          //System.out.println("Shifting " + (31-index));
          breaks |= 1<<index;
          count++;
        }
        if (count>=extend) break;
      }
      
      Sequence s = new Sequence(this.start, this.length, breaks);
      return s;
    }
    
    /*
    public int totalMass() {
      Iterator<Integer> it = iterator();
      int cumMass = 0;
      while (it.hasNext()) {
        cumMass += it.next();
      }
      return cumMass;
    }*/
    
    
    /**
     * Make sure all masses in this sequence are less than or equal than the limit
     * @param upperLimit the upper limit
     * @return true if all masses are under upper limit
     */
    /*
    public boolean checkMasses(int upperLimit) {
      Iterator<Integer> it = iterator();
      while (it.hasNext()) {
        if (it.next()>upperLimit) return false;
      }
      return true;
    }
    */
    
    @Override
    public int compareTo(Sequence other) {
      // lexicographical other
      Iterator<Integer> thisIt = iterator(), otherIt = other.iterator();
      while (true) {
        if (!thisIt.hasNext() && !otherIt.hasNext()) return 0;
        if (!thisIt.hasNext()) return -1; // other is bigger
        if (!otherIt.hasNext()) return 1; // this is bigger
        
        int thisMass = thisIt.next(), otherMass = otherIt.next();
        if (thisMass>otherMass) return 1;
        if (otherMass>thisMass) return -1;
      }
    }
    
    /**
     * Return the number of masses that 2 sequences have in common
     * @param other the other sequence
     * @return the number of masses that these 2 sequences have in common.
     */
    public int longestPrefix(Sequence other) {
      Iterator<Integer> thisIt = iterator(), otherIt = other.iterator();
      int count = 0;
      while (true) {
        if (!thisIt.hasNext() || !otherIt.hasNext()) return count;
        
        int thisMass = thisIt.next(), otherMass = otherIt.next();
        if (thisMass!=otherMass) return count;
        count++;
      }
    }

    @Override
    public Iterator<Integer> iterator() {
      return new MassIterator();
    }
    
    @Override
    public String toString() {
      StringBuffer sb = new StringBuffer();
      Iterator<Integer> it = iterator();
      while(it.hasNext()) {
        sb.append(it.next() + " ");
      }
      
      /*
      for (int i = 0; i < 5; i++) {
        sb.append((breaks>>>(31-i))&1);
      }*/
      return sb.toString();
    }
    
  }
  
  private ProteinFastaSequence sequence;
  
  public SequenceFactory(ProteinFastaSequence sequence) {
    this.sequence = sequence;
  }
  
  public Sequence getSequence(long start) {
    // limit the extension until we reach the first terminator
    int cumMass = 0;
    for (long i=start; i < 31; i++) {
      if (this.sequence.isTerminator(i)) {  
        return new Sequence(start, (int)(i-start), 0);
      }
      cumMass += this.sequence.getIntegerMass(start+i);
      if (cumMass > Constants.MAX_QUERY_MASS) break;
    }
    return new Sequence(start, 31, 0);
  }
  
  public Sequence getSequence(long start, int maxLength) {
    // limit the extension until we reach the first terminator
    int cumMass = 0;
    for (long i=start; i < maxLength; i++) {
      if (this.sequence.isTerminator(i)) {  
        return new Sequence(start, (int)(i-start), 0);
      }
      cumMass += this.sequence.getIntegerMass(start+i);
    }
    return new Sequence(start, maxLength, 0);
  }
  
  
  /**
   * 
   * @param sequences sorted sequences
   */
  private static void prune(ArrayList<Sequence> sequences) {
    //Collections.sort(sequences);
    Iterator<Sequence> it = sequences.iterator();
    Sequence current = it.next();
    while (it.hasNext()) {
      Sequence next = it.next();
      if (current.compareTo(next)==0) it.remove();
      current = next;
    }
  }
  
  
  public static void printRepetitionStats(ArrayList<Sequence> sequences, int start) {
    
    int i = start;
    
    //HashMap<Integer,Float> aveSizes = new HashMap<Integer,Float>();
    
    while (true) {
      
      Iterator<Sequence> it = sequences.iterator();
      Sequence prev = it.next();
      int mult = 1;
      boolean cont = false;
      // repetitiveness -> instances
      TreeMap<Integer,Integer> stats = new TreeMap<Integer,Integer>();
      while (it.hasNext()) {
        Sequence s = it.next();
        if (s.massCount >= i) {
          cont = true;
          int lcp = prev.longestPrefix(s);
          if (lcp==i) {
            mult++; // the is a k-mer
          }
          else {
            if (!stats.containsKey(mult)) {
              stats.put(mult, 0);
            }
            stats.put(mult, stats.get(mult)+1);
            
            prev = s; // reset the previous k-mer
            mult = 1;
          }
        }
        else {
          if (!stats.containsKey(mult)) {
            stats.put(mult, 0);
          }
          stats.put(mult, stats.get(mult)+1);
          
          prev = s;
          mult = 1;
        }
      }
      if (!stats.containsKey(mult)) {
        stats.put(mult, 0);
      }
      stats.put(mult, stats.get(mult)+1);
      
      // print our the stats for these k-mers
      System.out.print("k-mers of length " + i + " ");
      int cumCount = 0;
      int totalCount = 0;
      for (int key : stats.keySet()) {
        System.out.printf("%d:%d ", key, stats.get(key));
        if (key > 1) {
          cumCount += key * stats.get(key);
          totalCount += stats.get(key);
        }
      }
      System.out.println();
      //System.out.printf("%f\n", cumCount / (float)totalCount);
      
      i++;
      
      if (!cont) break;
    }
  }
  
  
  
  private static void expand(Sequence s, int index, ArrayList<Sequence> results) {
    int extend = 1;
    while (index+extend < s.massCount) {
      Sequence extra = s.derive(index, extend);
      if (extra.maxMass > Constants.MAX_GAP_MASS) break;
      results.add(extra);
      //System.out.println("  Added " + extra);
      extend++;
    }
  }
  
  /**
   * Generalized extension for i-unique sequences that allows for up to w
   * repeated sequences
   * @param sequences
   * @param index
   * @param w
   */
  public static void extendGap(ArrayList<Sequence> sequences, int index, int w) {
    ArrayList<Sequence> extras = new ArrayList<Sequence>();
    
    Iterator<Sequence> it = sequences.iterator();
    ArrayList<Sequence> chunk = new ArrayList<Sequence>();
    ArrayList<Integer> lcps = new ArrayList<Integer>();
    chunk.add(it.next());
    
    
    while (it.hasNext()) {
      
      Sequence next = it.next();
      
      int lcp = chunk.get(chunk.size()-1).longestPrefix(next);
      
      //if (lcp-1 == index) {
      if (lcp-1 >= index) {
        chunk.add(next);
        lcps.add(lcp);
      }
      else {
        // expand all sequences, reset chunk
        if (chunk.size() > w) {
          boolean expand = false;
          for (int commonPrefix : lcps) {
            if (commonPrefix - 1 <= index) expand = true; 
          }
          if (expand || index==0) { // always extend at 0
            //System.out.println("Expanding " + index + " chunk size " + chunk.size());
            for (Sequence subS : chunk) {
              //System.out.println(subS);
              expand(subS, index, extras);
            }
          }
          else {
            /*
            System.out.println("Skipping " + index + " chunk size " + chunk.size());
            for (Sequence subS : chunk) {
              System.out.println(subS);
            }
            System.out.println();
            */
          }
        }
        else {
          //System.out.println("Skipping " + index + " w-chunk size " + chunk.size());
          //for (Sequence subS : chunk) {
          //  System.out.println(subS);
          //}
        }
        chunk.clear();
        lcps.clear();
        chunk.add(next);
      }
      
    }
    
    if (chunk.size() > w) {
      boolean expand = false;
      for (int commonPrefix : lcps) {
        if (commonPrefix - 1 <= index) expand = true; 
      }
      if (expand || index==0) { 
        for (Sequence subS : chunk) {
          expand(subS, index, extras);
        }
      }
    }
    sequences.addAll(extras);
    Collections.sort(sequences);
    prune(sequences);
  }

  /**
   * Generalized extension for i-unique sequences that allows for up to w.
   * NO - not optimized
   * repeated sequences
   * @param sequences
   * @param index
   * @param w
   */
  public static void extendGapNO(ArrayList<Sequence> sequences, int index, int w) {
    ArrayList<Sequence> extras = new ArrayList<Sequence>();
    
    Iterator<Sequence> it = sequences.iterator();
    ArrayList<Sequence> chunk = new ArrayList<Sequence>();
    chunk.add(it.next());
    
    while (it.hasNext()) {
      
      Sequence next = it.next();
      
      int lcp = chunk.get(chunk.size()-1).longestPrefix(next);
      
      //if (lcp-1 == index) {
      if (lcp-1 >= index) {
        chunk.add(next);
      }
      else {
        // expand all sequences, reset chunk
        if (chunk.size() > w) {
          //System.out.println("Expanding " + index + " chunk size " + chunk.size());
          for (Sequence subS : chunk) {
            //System.out.println(subS);
            expand(subS, index, extras);
          }
        }
        else {
          //System.out.println("Skipping " + index + " w-chunk size " + chunk.size());
          //for (Sequence subS : chunk) {
          //  System.out.println(subS);
          //}
        }
        chunk.clear();
        chunk.add(next);
      }
      
    }
    
    if (chunk.size() > w) {
      for (Sequence subS : chunk) {
        expand(subS, index, extras);
      }
    }
    sequences.addAll(extras);
    Collections.sort(sequences);
    prune(sequences);
  }


  public static void main(String[] args) {
    String userHome = System.getProperty("user.home");
    String fastaFile;
    //String outFile;
    
    int K = 15;
    int W = 1;
    
    fastaFile = userHome+"/Data/Databases/test.fasta";
    fastaFile = userHome+"/Data/Databases/tiny.fasta";
    //fastaFile = userHome+"/Data/Databases/small.fasta";
    //fastaFile = userHome+"/Data/Databases/random1000.fasta";
    fastaFile = userHome+"/Data/Databases/random100000.fasta";
    fastaFile = userHome+"/Data/Databases/random200000.fasta";
    //fastaFile = userHome+"/Data/Databases/random500000.fasta";
    //fastaFile = userHome+"/Data/Databases/random1000000.fasta";
    //fastaFile = userHome+"/Data/Databases/half.fasta";
    //fastaFile = userHome+"/Data/Databases/million.fasta";
    //fastaFile = userHome+"/Data/Databases/debug.fasta";  
    fastaFile = userHome+"/Data/Databases/ShewDB/SOne_uniprot_plus_contaminants.fasta";
    //fastaFile = userHome+"/Data/Databases/Smuelleri.fasta";
    
    ProteinFastaSequence sequence = new ProteinFastaSequence(fastaFile);
    SequenceFactory sf = new SequenceFactory(sequence);
    
    ArrayList<Sequence> sequences = new ArrayList<Sequence>();
    // create all K-mers
    //for (long i=0; i<sequence.getSize(); i++) {
    int MAX_SEQUENCE_COUNT = 300000;
    for (long i=0; i<MAX_SEQUENCE_COUNT; i++) {
      if (sequence.isTerminator(i)) continue;
      //Sequence s = sf.getSequence(i);
      Sequence s = sf.getSequence(i, K);
      sequences.add(s);
    }
    
    Collections.sort(sequences);
    prune(sequences);
    System.out.println("Unique sequence count " + sequences.size());

    //printRepetitionStats(sequences, 1);
    
    //Collections.sort(sequences);
    for (int i = 0; i<K; i++) {
      System.out.println(i + " " + sequences.size());
      extendGapNO(sequences, i, W);
      //extendGap(sequences, i, W);
      //break;
    }
    
    System.out.println("Total sequences " + sequences.size());
    
    // build a keyword tree out of it and collect stats
    ArrayList<ArrayList<Integer>> queries = new ArrayList<ArrayList<Integer>>();
    for (Sequence s : sequences) {
      queries.add(s.toArray());
    }
    KeywordTree kt = new KeywordTree(queries, false);
    System.out.println(kt.collectStats().toString());
    
  }

}
