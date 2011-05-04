package suffixgraph.test;

import java.util.ArrayList;
import java.util.Random;

import msutil.GappedPeptide;

import sequences.FastaSequence;
import suffixgraph.graphs.IntegerGappedSuffixGraph;


/**
 * Debugging and testing code for the GappedSuffixGraph class.
 * @author jung
 *
 */
public class IntegerGappedSuffixGraphTest {
  
  
  private static String generateRandomGap(Random gen, String peptide) {
    String finalPeptide = "";
    for (int i=0; i<peptide.length();) {
      int gapLen = gen.nextInt(3) + 1;
      finalPeptide += "[" + peptide.substring(i, Math.min(i+gapLen, peptide.length())) + "]";
      i += gapLen;
    }
    return finalPeptide;    
  }
  
  
  /**
   * Exhaustively test all sequences in the graph.
   */
  private static void allTest() {
    String fastaFile;
    String userHome = System.getProperty("user.home");
    
    fastaFile = userHome+"/Data/Databases/test.fasta";
    fastaFile = userHome+"/Data/Databases/small.fasta";
    fastaFile = userHome+"/Data/Databases/medium.fasta";
    //fastaFile = userHome+"/Data/Databases/large.fasta";
    //fastaFile = userHome+"/Data/Databases/uniprot_sprot.fasta";
    //fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    
    String[] tokens = fastaFile.split("/");
    String graphFile = userHome+"/Desktop/"+tokens[tokens.length-1].replaceAll(".fasta$", "")+".igsg";
    
    long time = System.currentTimeMillis();
    FastaSequence sequence = new FastaSequence(fastaFile);
    System.out.println("-- Loading fasta file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    System.out.println("--- Number of characters in fasta file: " + sequence.getSize());
    
    time = System.currentTimeMillis();
    IntegerGappedSuffixGraph gsg = new IntegerGappedSuffixGraph(sequence, graphFile);
    System.out.println("-- Loading SuffixGraph file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    long usedMem = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024);
    System.out.println("---- Memory used for the gapped graph " + usedMem + "MB");
    
    int queryCount = 0;
    int cumTime = 0;
    Random gen = new Random();
    // query all possible sequences of a certain length
    for (int i=0; i < sequence.getSize(); i++) {
      String queryStr = sequence.getSubsequence(i, i+15);
      if (queryStr.contains("_")) {
        continue;
      }
      
      try {
        GappedPeptide query = new GappedPeptide(queryStr);
        time = System.currentTimeMillis();
        ArrayList<String> matches = gsg.search(query, query.getCount());
        cumTime += System.currentTimeMillis() - time;
        queryCount++;
        if (matches == null) {
          System.err.println("Not found " + queryStr);
        }

        String gappedQuery = generateRandomGap(gen, queryStr); 
        //System.out.println(gappedQuery);
        // query a random gapped version of this peptide
        matches = gsg.search(new GappedPeptide(gappedQuery));
        time = System.currentTimeMillis();
        cumTime += System.currentTimeMillis() - time;
        queryCount++;
        if (matches == null) {
          System.err.println("Not found " + gappedQuery);
        }
        
      }
      catch (NullPointerException e) {
      }
     
    }
    System.out.println("-- Searching " + queryCount + " queries within " + cumTime/1000.0 + "s");
    
  }
  
  
  public static void main(String[] args) {
    allTest();
  }
  
}
