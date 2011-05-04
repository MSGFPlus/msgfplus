package suffixgraph.test;

import java.util.ArrayList;
import java.util.Random;

import msutil.GappedPeptide;

import sequences.FastaSequence;
import suffixgraph.graphs.CompositionGappedSuffixGraph;


/**
 * Debugging and testing code for the GappedSuffixGraph class.
 * @author jung
 *
 */
public class CompositionGappedSuffixGraphTest {
  
  
  private static String generateRandomGap(Random gen, String peptide) {
    String finalPeptide = "";
    for (int i=0; i<peptide.length();) {
      int gapLen = gen.nextInt(2) + 1;
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
    
    //fastaFile = userHome+"/Data/Databases/test.fasta";
    //fastaFile = userHome+"/Data/Databases/small.fasta";
    fastaFile = userHome+"/Data/Databases/medium.fasta";
    //fastaFile = userHome+"/Data/Databases/million.fasta";
    //fastaFile = userHome+"/Data/Databases/large.fasta";
    //fastaFile = userHome+"/Data/Databases/uniprot_sprot.fasta";
    //fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    //fastaFile = userHome+"/Data/Databases/ShewDB/SOne_proteins_withContams.fasta";
    
    String graphFile = fastaFile.replaceAll(".fasta$", "")+".cgsg";
    
    long time = System.currentTimeMillis();
    FastaSequence sequence = new FastaSequence(fastaFile);
    System.out.println("-- Loading fasta file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    System.out.println("--- Number of characters in fasta file: " + sequence.getSize());
    
    time = System.currentTimeMillis();
    CompositionGappedSuffixGraph gsg = new CompositionGappedSuffixGraph(sequence, graphFile);
    System.out.println("-- Loading SuffixGraph file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    System.out.println(gsg);
    
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
