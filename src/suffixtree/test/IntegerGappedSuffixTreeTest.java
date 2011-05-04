package suffixtree.test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;

import sequences.FastaSequence;
import sequences.ProteinFastaSequence;
import suffixtree.Constants;
import suffixtree.matches.ExactMatchObject;
import suffixtree.trees.deprecated.IntegerGappedSuffixTree;

public class IntegerGappedSuffixTreeTest {
 
  public static String gappedQueryString(FastaSequence sequence, ArrayList<Integer> query) {
    StringBuffer sb = new StringBuffer();
    sb.append("{");
    for(int qe : query) {
      sb.append(qe + ", ");
    }
    return sb.toString();
  }
  
  
  public static long query(IntegerGappedSuffixTree st, ArrayList<Integer> query) {
    long time = System.currentTimeMillis();
    //System.err.println("Querying: " + gappedQueryString(st.getSequence(), query));
    HashSet<ExactMatchObject> matches = new HashSet<ExactMatchObject>();
    st.search(query, matches);
    time = System.currentTimeMillis() - time;
    
    if (matches.size()==0) {
      System.err.print("Querying: " + gappedQueryString(st.getSequence(), query));
      System.err.println(" - Not found!");
      System.exit(-1);
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
  
  
  private static void queryAll(ProteinFastaSequence sequence, IntegerGappedSuffixTree st) {
    Random r = new Random();
    
    int queryCount = 0;
    long cumTime = 0;
    // query all the items
    for (int start=0; start < sequence.getSize(); start++) {
      if (sequence.isTerminator(start)) continue;
      
      // for each start position, find the end
      int end = start;
      int cumMass = sequence.getIntegerMass(start);
      for (int i=start+1; i < sequence.getSize(); i++) {
         if (sequence.isTerminator(i) || cumMass>Constants.MAX_QUERY_MASS) {
          end = i;
          break;
        }
        cumMass += sequence.getIntegerMass(i);
      }
      
      if (cumMass<Constants.MIN_QUERY_MASS) continue;
      
      //System.err.println("-- Querying... " + sequence.toString(start, end));
      
      // intact query
      ArrayList<Integer> query = new ArrayList<Integer>();
      for (int i=start; i<end; i++) {
        query.add(sequence.getIntegerMass(i));
      }
      queryCount++;
      cumTime += query(st, query);
      
      //gapped query
      ArrayList<Integer> gappedQuery = new ArrayList<Integer>();
      for (int i=start; i<end;) {
        int gapSize = r.nextInt(Constants.MAX_GAP)+1;
        if (i+gapSize >= end) break;
        
        cumMass = 0;
        int j;
        for (j=0; j<gapSize && cumMass+186<=Constants.MAX_GAP_MASS; j++) 
          cumMass += sequence.getIntegerMass(i+j); 
        gappedQuery.add(cumMass);
        i += j;
      }        
      
      if (gappedQuery.size()==0) continue;
      
      queryCount++;
      cumTime += query(st, gappedQuery);
    }
    
    System.out.printf("-- %d queries in %.2f seconds\n", queryCount, cumTime/1000.0);
    System.out.printf("-- Average %.2f ms per 1000 query", 1000.0*cumTime/queryCount);
  }
  
  public static void main(String[] args) {
    String fastaFile;
    String userHome = System.getProperty("user.home");
    
    fastaFile = userHome+"/Data/Databases/test.fasta";
    fastaFile = userHome+"/Data/Databases/tiny.fasta";
    fastaFile = userHome+"/Data/Databases/small.fasta";
    //fastaFile = userHome+"/Data/Databases/medium.fasta";
    //fastaFile = userHome+"/Data/Databases/repeat.fasta";
    //fastaFile = userHome+"/Data/Databases/half.fasta";
    //fastaFile = userHome+"/Data/Databases/million.fasta";
    //fastaFile = userHome+"/Data/Databases/large.fasta";
    //fastaFile = userHome+"/Data/Databases/uniprot_sprot.fasta";
    //fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    //fastaFile = userHome+"/Data/Databases/ecoli.fasta";
    //fastaFile = userHome+"/Data/Databases/ShewDB/SOne_uniprot_plus_contaminants.fasta";
    //fastaFile = userHome+"/Data/Databases/random100000.fasta";
    //fastaFile = userHome+"/Data/Databases/random500000.fasta";
    //fastaFile = userHome+"/Data/Databases/random1000000.fasta";
    //fastaFile = userHome+"/Data/Databases/Mgenitalium.fasta";
    
    //String graphFile = fastaFile.replaceAll(".fasta$", "")+".st";
   
    long time = System.currentTimeMillis();
    ProteinFastaSequence sequence = new ProteinFastaSequence(fastaFile, sequences.Constants.AMINO_ACIDS_18);
    System.out.println("-- Loading fasta file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    System.out.println("--- Number of characters in fasta file: " + sequence.getSize());
    
    time = System.currentTimeMillis();
    IntegerGappedSuffixTree st = new IntegerGappedSuffixTree(sequence);
    System.out.println("-- Loading SuffixTree file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    System.out.println(st.collectStats());
    //System.out.print(st);

    queryAll(sequence, st);
    
    /*
    //Querying... ILGADELVMSPIPTTDVQPKVTFDINSEVSSGPLYLNPVEMAGVKYLQL
    int[] qEdgesInts = {113, 113, 128, 115, 129, 113};
    ArrayList<QueryEdge> query = new ArrayList<QueryEdge>();
    for (int i : qEdgesInts) {
      query.add(new QueryEdge(i));
    }
    
    HashSet<Integer> matches = new HashSet<Integer>(); 
    st.search(query, matches);
    if (matches.size()>0) {
      System.err.println("Matches [" + matches.size() + "] at: ");
      for (int i : matches) {
        System.err.print(i + " ");
      }
      System.err.println();
    }
    else {
      System.err.println("Not found");
    }
    */
  }
}
