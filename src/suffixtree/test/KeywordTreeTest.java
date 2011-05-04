package suffixtree.test;

import java.util.ArrayList;

import sequences.ProteinFastaSequence;
//import sequences.ProteinFastaSequences;
import suffixtree.matches.ExactMatchObject;
import suffixtree.trees.KeywordTree;



/**
 * Benchmark and validator class of the KeywordTree
 * @author jung
 *
 */
public class KeywordTreeTest {
  
  
  public static void queryAll(ProteinFastaSequence sequence, KeywordTree kt, int maxItems, String outfile) {
    if (outfile==null) {
      ArrayList<ExactMatchObject> matches = new ArrayList<ExactMatchObject>();
      kt.match(sequence, matches);
    }
    else {
      //kt.exactMatch(sequence, outfile);
    }
  }

  
  
  public static void main(String[] args) {
    String userHome = System.getProperty("user.home");
    String fastaFile;
    String outFile;
    
    int maxItems = 100000;
    
    fastaFile = userHome+"/Data/Databases/test.fasta";
    fastaFile = userHome+"/Data/Databases/small.fasta";
    fastaFile = userHome+"/Data/Databases/medium.fasta";
    //fastaFile = userHome+"/Data/Databases/half.fasta";
    fastaFile = userHome+"/Data/Databases/million.fasta";
    //fastaFile = userHome+"/Data/Databases/large.fasta";
    //fastaFile = userHome+"/Data/Databases/uniprot_sprot.fasta";
    //fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    //fastaFile = userHome+"/Data/Databases/ShewDB/SOne_proteins_withContams.fasta";
    //fastaFile = userHome+"/Data/Databases/miniDir";
    //fastaFile = userHome+"/Data/Databases/Human.54";
    
    outFile = userHome+"/Desktop/results.txt";
    //outFile = null;
    
    long time = System.currentTimeMillis();
    ProteinFastaSequence sequence = new ProteinFastaSequence(fastaFile, sequences.Constants.AMINO_ACIDS_18);
    //ProteinFastaSequence sequence = new ProteinFastaSequences(fastaFile, FastaSequence.AminoAcids18, false);
    System.out.println("-- Loading fasta file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    System.out.println("-- Number of characters in fasta file: " + sequence.getSize());
    
    time = System.currentTimeMillis();
    ArrayList<ArrayList<Integer>> queries = new ArrayList<ArrayList<Integer>>();
    SequenceGenerator.generateRandomCorrectQueries(sequence, queries, maxItems);
    KeywordTree kt = new KeywordTree(queries);
    System.out.println("-- Done building the keyword tree in " + (System.currentTimeMillis() - time)/1000.0 + "s");
    
    time = System.currentTimeMillis();
    queryAll(sequence, kt, maxItems, outFile);
    
    System.out.println("-- Done making queries in " + (System.currentTimeMillis() - time)/1000.0 + "s");
  }
}
