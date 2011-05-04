package suffixtree.test;

import java.util.ArrayList;
import java.util.HashSet;

import sequences.ProteinFastaSequence;
import suffixtree.matches.ExactMatchObject;
import suffixtree.trees.FailingKeywordTree;



/**
 * Benchmark and validator class of the KeywordTree
 * @author jung
 *
 */
public class FailingKeywordTreeTest {
  
  public static void queryAll(ProteinFastaSequence sequence, FailingKeywordTree fkt, int maxItems) {
    HashSet<ExactMatchObject> matches = new HashSet<ExactMatchObject>();
    long time = System.nanoTime();
    fkt.feed(sequence, matches);
    time = System.nanoTime() - time;
    double aveTime = time / (1000.0*maxItems);
    System.out.printf("Querying %.2f ms per 1000 queries\n", aveTime);
    
    System.out.println("Total matches: " + matches.size() + ". Total queries: " + maxItems);
  }

  
  
  public static void main(String[] args) {
    String userHome = System.getProperty("user.home");
    String fastaFile;
    
    int maxItems = 1000000;
    
    fastaFile = userHome+"/Data/Databases/test.fasta";
    fastaFile = userHome+"/Data/Databases/debug.fasta";
    fastaFile = userHome+"/Data/Databases/small.fasta";
    fastaFile = userHome+"/Data/Databases/medium.fasta";
    //fastaFile = userHome+"/Data/Databases/half.fasta";
    //fastaFile = userHome+"/Data/Databases/million.fasta";
    //fastaFile = userHome+"/Data/Databases/large.fasta";
    //fastaFile = userHome+"/Data/Databases/uniprot_sprot.fasta";
    //fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    //fastaFile = userHome+"/Data/Databases/ShewDB/SOne_proteins_withContams.fasta";
    
    long time = System.currentTimeMillis();
    ProteinFastaSequence sequence = new ProteinFastaSequence(fastaFile, sequences.Constants.AMINO_ACIDS_18);
    System.out.println("-- Loading fasta file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    System.out.println("-- Number of characters in fasta file: " + sequence.getSize());
    
    time = System.currentTimeMillis();
    ArrayList<ArrayList<Integer>> queries = new ArrayList<ArrayList<Integer>>();
    SequenceGenerator.generateRandomCorrectQueries(sequence, queries, maxItems);
    FailingKeywordTree fkt = new FailingKeywordTree(queries);
    System.out.println("-- Done building the keyword tree in " + (System.currentTimeMillis() - time)/1000.0 + "s");
    
    time = System.currentTimeMillis();
    queryAll(sequence, fkt, maxItems);
    System.out.println("-- Done making queries in " + (System.currentTimeMillis() - time)/1000.0 + "s");
  }
}
