package suffixtree.test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;

import sequences.FastaSequence;
import suffixtree.Constants;
import suffixtree.edges.ByteEdge;
import suffixtree.trees.HashedGappedSuffixTree;

public class HashedGappedSuffixTreeTest {
  
  public static String gappedQueryString(FastaSequence sequence, ArrayList<ByteEdge> query) {
    StringBuffer sb = new StringBuffer();
    for(ByteEdge qe : query) {
      sb.append("[" + sequence.toString(qe.toByteArray()) + "," + qe.getLabel() + "] ");
    }
    return sb.toString();
  }
  
  public static long query(HashedGappedSuffixTree st, ArrayList<ByteEdge> query) {
    long time = System.currentTimeMillis();
    //System.out.println("Querying: " + gappedQueryString(st.getSequence(), query));
    HashSet<Integer> matches = new HashSet<Integer>();
    st.search(query, matches);
    time = System.currentTimeMillis() - time;
    
    if (matches.size()==0) {
      System.out.print("Querying: " + gappedQueryString(st.getSequence(), query));
      System.out.println(" - Not found!");
      return -1;
    }
    else {
      /*
      if (matches.size() > 10 || time > 1000) {
        System.out.println(" in " + time + " milisecs");
        System.out.print(" - [" + matches.size() + "] positions: ");
        for (int i : matches) {
          System.out.println("     " + i + " : " + st.getSequence().toString(i, i+MIN_QUERY_LENGTH));
        }
        System.out.println();
      }*/
    }
    
    return time;
  }
  
  public static void queryAll(FastaSequence sequence, HashedGappedSuffixTree st) {
    Random r = new Random();
    
    int queryCount = 0;
    long cumTime = 0;
    // query all the items
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
      
      end = Math.min(end, start+Constants.MAX_QUERY_CHAR);
      if (end-start >= Constants.MIN_QUERY_CHAR) {
      
        // intact query
        ArrayList<ByteEdge> query = new ArrayList<ByteEdge>();
        for (int i=start; i<end; i++) {
          query.add(new ByteEdge(sequence.getByteAt(i)));
        }
        queryCount++;
        cumTime += query(st, query);
        
        //gapped query
        ArrayList<ByteEdge> gappedQuery = new ArrayList<ByteEdge>();
        for (int i=start; i<end;) {
          int gapSize = r.nextInt(Constants.MAX_GAP)+1;
          if (i+gapSize >= end) break;
          gappedQuery.add(new ByteEdge(sequence.getBytes(i, i+gapSize)));
          i += gapSize;
        }        
        
        if (gappedQuery.size()==0) continue;
        
        queryCount++;
        cumTime += query(st, gappedQuery);
      }
      
      if (start % 100000 == 0) {
        System.out.println("---- Queried approximately " + start + " items.");
      }
    }
    
    System.out.printf("-- %d queries in %.2f seconds\n", queryCount, cumTime/1000.0);
    System.out.printf("-- Average %.2f ms per 1000 query", 1000.0*cumTime/queryCount);
  }
  
  
  public static void main(String[] args) {
    String userHome = System.getProperty("user.home");
    String fastaFile;
    
    fastaFile = userHome+"/Data/Databases/test.fasta";
    fastaFile = userHome+"/Data/Databases/small.fasta";
    fastaFile = userHome+"/Data/Databases/medium.fasta";
    //fastaFile = userHome+"/Data/Databases/half.fasta";
    //fastaFile = userHome+"/Data/Databases/million.fasta";
    //fastaFile = userHome+"/Data/Databases/large.fasta";
    //fastaFile = userHome+"/Data/Databases/uniprot_sprot.fasta";
    //fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    fastaFile = userHome+"/Data/Databases/ShewDB/SOne_uniprot_plus_contaminants.fasta";
    
    long time = System.currentTimeMillis();
    FastaSequence sequence = new FastaSequence(fastaFile);
    System.out.println("-- Loading fasta file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    System.out.println("--- Number of characters in fasta file: " + sequence.getSize());
    
    time = System.currentTimeMillis();
    HashedGappedSuffixTree hgst = new HashedGappedSuffixTree(sequence);
    System.out.println("--- Done building the tree in " + (System.currentTimeMillis() - time)/1000.0 + "s");
    queryAll(sequence, hgst);
  }
  

}
