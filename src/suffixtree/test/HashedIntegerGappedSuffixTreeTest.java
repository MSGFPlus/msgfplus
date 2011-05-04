package suffixtree.test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;

import sequences.FastaSequence;
import sequences.ProteinFastaSequence;
import suffixtree.Constants;
import suffixtree.matches.ExactMatchObject;
import suffixtree.trees.HashedIntegerGappedSuffixTree;

public class HashedIntegerGappedSuffixTreeTest {

  public static String gappedQueryString(FastaSequence sequence, ArrayList<Integer> query) {
    StringBuffer sb = new StringBuffer();
    for(int qe : query) {
      sb.append(qe + ", ");
    }
    return sb.toString();
  }
  
  public static long query(HashedIntegerGappedSuffixTree st, ArrayList<Integer> query) {
    long time = System.nanoTime();
    //System.out.println("Querying: " + gappedQueryString(st.getSequence(), query));
    HashSet<ExactMatchObject> matches = new HashSet<ExactMatchObject>();
    st.search(query, matches);
    time = System.nanoTime() - time;
    
    if (matches.size()==0) {
      System.out.print("Querying: " + gappedQueryString(st.getSequence(), query));
      System.out.println(" - Not found!");
      //System.exit(-1);
      return -1;
    }
    else {
      /*
      System.out.println(" in " + time / 1000000.0 + " milisecs");
      System.out.print(" - [" + matches.size() + "] positions: ");
      for (int i : matches) {
        System.out.println("     " + i + " : " + st.getSequence().toString(i, i+Constants.MIN_QUERY_LENGTH));
      }
      System.out.println();
      */
    }
    
    return time;
  }
  
  public static void queryAll(ProteinFastaSequence sequence, HashedIntegerGappedSuffixTree st) {
    Random r = new Random();
    
    int queryCount = 0;
    long cumTime = 0;
    // query all the items
    for (int start=1; start < sequence.getSize(); start++) {
      //System.err.println("-- Querying... " + sequence.toString(start, end));
        
      // intact query
      ArrayList<Integer> masses = SequenceGenerator.generateIntegerMassSet(sequence, start, Constants.MIN_QUERY_MASS, Constants.MAX_QUERY_MASS);
      if (masses != null) {
        queryCount++;
        long time = query(st, masses);
        if (time < 0) {
          System.out.println("-- Not found... " + sequence.getSubsequence(start, 10+start));
        }
        cumTime += time;
      }
       
      //gapped query
      masses = SequenceGenerator.generateRandomIntegerMassSet(sequence, start, Constants.MIN_QUERY_MASS, Constants.MAX_QUERY_MASS, Constants.MAX_GAP_MASS, r);
      if (masses != null) {
        queryCount++;
        long time = query(st, masses);
        if (time < 0) {
          System.out.println("-- Not found... " + sequence.getSubsequence(start, 10+start));
        }
        cumTime += time;
      }
      
      if (start % 100000 == 0) {
        System.out.println("---- Queried approximately " + 2*start + " items.");
      }
    }
    
    System.out.printf("-- %d queries in %.2f seconds\n", queryCount, cumTime/1000000000.0);
    System.out.printf("-- Average %.2f ms per 1000 query", cumTime/(1000.0*queryCount));
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
    fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    //fastaFile = userHome+"/Data/Databases/ShewDB/SOne_uniprot_plus_contaminants.fasta";
    //fastaFile = userHome+"/Data/Databases/ShewDB/SOne_uniprot_plus_contaminants_ALL.fasta";
    
    long time = System.currentTimeMillis();
    ProteinFastaSequence sequence = new ProteinFastaSequence(fastaFile, sequences.Constants.AMINO_ACIDS_18);
    System.out.println("-- Loading fasta file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    System.out.println("--- Number of characters in fasta file: " + sequence.getSize());
    
    time = System.currentTimeMillis();
    HashedIntegerGappedSuffixTree higst = new HashedIntegerGappedSuffixTree(sequence);
    System.out.println("--- Done building the tree in " + (System.currentTimeMillis() - time)/1000.0 + "s");
    queryAll(sequence, higst);
  
    
    // EINVTPGAISHKISTIEDFIGKKVFERGSRRVTL
    //int[] masses = {129, 113, 114, 99, 101, 97, 128, 200, 378, 87, 101, 242, 262, 298, 128, 99, 276, 300, 411, 214, 424, 432, 406, 313, 226, 128};
    //int[] masses = {129, 113, 411, 128, 200, 378, 87, 101, 242, 262, 298, 128, 99, 276, 300, 411, 214, 424, 432, 406, 313, 226, 128};
    /*
    int[] masses = {129, 113, 411, 128, 113+87, 137, 128, 113, 87, 101, 242, 262, 298, 128, 99, 276, 300, 411, 214, 424, 432, 406, 313, 226, 128};
    ArrayList<QueryEdge> query = new ArrayList<QueryEdge>();
    for (int mass : masses) {
      query.add(new QueryEdge(mass));
    }
    query(higst, query);
    */
  }
  

}
