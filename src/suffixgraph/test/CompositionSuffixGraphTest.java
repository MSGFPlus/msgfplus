package suffixgraph.test;

import java.util.ArrayList;

import msutil.GappedPeptide;
import sequences.FastaSequence;
import suffixgraph.graphs.CompositionSuffixGraph;

/**
 * Debugging and testing code for the GappedSuffixGraph class.
 * @author jung
 *
 */
public class CompositionSuffixGraphTest {

  /**
   * Tester method.
   */
  private static void testAll() {
    String fastaFile;
    String userHome = System.getProperty("user.home");
    
    fastaFile = userHome+"/Data/Databases/test.fasta";
    //fastaFile = userHome+"/Data/Databases/small.fasta";
    //fastaFile = userHome+"/Data/Databases/medium.fasta";
    //fastaFile = userHome+"/Data/Databases/large.fasta";
    //fastaFile = userHome+"/Data/Databases/uniprot_sprot.fasta";
    //fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    
    String[] tokens = fastaFile.split("/");
    String graphFile = userHome+"/Desktop/"+tokens[tokens.length-1].replaceAll(".fasta$", "")+".csg";
    
    long time = System.currentTimeMillis();
    FastaSequence sequence = new FastaSequence(fastaFile);
    System.out.println("-- Loading fasta file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    System.out.println("--- Number of characters in fasta file: " + sequence.getSize());
    
    time = System.currentTimeMillis();
    CompositionSuffixGraph st = new CompositionSuffixGraph(sequence, graphFile);
    System.out.println("-- Loading SuffixGraph file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    System.out.println(st);
    
    time = System.currentTimeMillis();
    long cumTime = 0;
    int queryCount = 0;
    for (int i=0; i < sequence.getSize(); i++) {
      String queryStr = sequence.getSubsequence(i, i+15);
      if (queryStr.contains("_")) {
        continue;
      }
      
      try {
        GappedPeptide query = new GappedPeptide(queryStr);
        time = System.currentTimeMillis();
        ArrayList<String> matches = st.search(query.getCompositions());
        cumTime += System.currentTimeMillis() - time;
        queryCount++;
        if (matches == null) {
          System.err.println("Not found " + queryStr);
        }
      }
      catch (NullPointerException e) {
      }
     
    }
    System.out.println("-- Searching " + queryCount + " queries within " + cumTime/1000.0 + "s");
  }
  
  
  public static void main(String[] args) {
    testAll();
  }
  
}