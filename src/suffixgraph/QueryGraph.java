package suffixgraph;

import java.util.ArrayList;
import java.util.Scanner;

import msutil.GappedPeptide;

import sequences.FastaSequence;
import suffixgraph.graphs.IntegerGappedSuffixGraph;

/**
 * Interactive querying of a database
 * @author jung
 *
 */
public class QueryGraph {
  
  /**
   * Exhaustively test all sequences in the graph.
   */
  public static void main(String[] args) {
    String fastaFile;
    String userHome = System.getProperty("user.home");
    
    String graphFile = userHome+"/Desktop/generic.graph";
    
    fastaFile = userHome+"/Data/Databases/test.fasta";
    fastaFile = userHome+"/Data/Databases/small.fasta";
    //fastaFile = userHome+"/Data/Databases/medium.fasta";
    //fastaFile = userHome+"/Data/Databases/large.fasta";
    //fastaFile = userHome+"/Data/Databases/uniprot_sprot.fasta";
    //fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    
    long time = System.currentTimeMillis();
    FastaSequence sequence = new FastaSequence(fastaFile);
    System.out.println("-- Loading fasta file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    System.out.println("--- Number of characters in fasta file: " + sequence.getSize());
    
    time = System.currentTimeMillis();
    IntegerGappedSuffixGraph gsg = new IntegerGappedSuffixGraph(sequence, graphFile);
    System.out.println("-- Loading SuffixGraph file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    long usedMem = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024);
    System.out.println("---- Memory used for the gapped graph " + usedMem + "MB");

    String queryLine = null;
    do {
      System.out.print("Enter query: ");
      Scanner in = new Scanner(System.in);
      queryLine = in.nextLine().trim();
      ArrayList<String> matches = gsg.search(new GappedPeptide(queryLine));
      if (matches != null) {
        for (String match : matches) {
          System.out.println(match + " found");
        }
      }
      else {
        System.out.println(queryLine + " not found");
      }
      
    }
    while (true);
  }
}
