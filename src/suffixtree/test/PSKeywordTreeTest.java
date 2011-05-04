package suffixtree.test;

import java.util.ArrayList;

import sequences.MassSequence;
import sequences.ProteinFastaSequence;
//import sequences.ProteinFastaSequences;
import suffixtree.matches.MatchObject;
import suffixtree.matches.MegaPrefixSuffixMatchObject;
import suffixtree.trees.ComplexKeywordTree;
import suffixtree.trees.PSKeywordTree;



/**
 * Benchmark and validator class of the KeywordTree
 * @author jung
 *
 */
public class PSKeywordTreeTest {
  
  
  public static void queryAll(ProteinFastaSequence sequence, ComplexKeywordTree kt, int maxItems, String outfile) {
    if (outfile==null) {
      //ArrayList<MatchObject> matches = new ArrayList<MatchObject>();
      //kt.feed(sequence, matches);
    }
    else {
      //kt.populate(sequence, outfile);
      //kt.reversePopulate(sequence, outfile);
    }
  }

  
  public static void main(String[] args) {
    String userHome = System.getProperty("user.home");
    String fastaFile;
    //String outFile;
    
    int maxItems = 0;
    
    fastaFile = userHome+"/Data/Databases/test.fasta";
    fastaFile = userHome+"/Data/Databases/small.fasta";
    fastaFile = userHome+"/Data/Databases/debug.fasta";
    //fastaFile = userHome+"/Data/Databases/single.fasta";
    //fastaFile = userHome+"/Data/Databases/medium.fasta";
    //fastaFile = userHome+"/Data/Databases/half.fasta";
    //fastaFile = userHome+"/Data/Databases/million.fasta";
    //fastaFile = userHome+"/Data/Databases/large.fasta";
    //fastaFile = userHome+"/Data/Databases/uniprot_sprot.fasta";
    //fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
    //fastaFile = userHome+"/Data/Databases/ShewDB/SOne_uniprot_plus_contaminants.fasta";
    //fastaFile = userHome+"/Data/Databases/miniDir";
    //fastaFile = userHome+"/Data/Databases/Human.54";
    
    //outFile = userHome+"/Desktop/results.txt";
    //outFile = null;
    
    long time = System.currentTimeMillis();
    MassSequence sequence = new ProteinFastaSequence(fastaFile, sequences.Constants.AMINO_ACIDS_18);
    //ProteinFastaSequence sequence = new ProteinFastaSequences(fastaFile, FastaSequence.AminoAcids18, false);
    System.out.println("-- Loading fasta file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    System.out.println("-- Number of characters in fasta file: " + sequence.getSize());
    
    /*
    for (int index = 0; index < sequence.getSize(); index++) {
      System.out.printf("%c %d\n", sequence.getCharAt(index), sequence.getIntegerMass(index));  
    }
    */
    
    
    // this query is in the debug sequence
    //int[] query = {131, 128, 147, 113, 99, 113, 113, 147, 114, 113, 113, 103, 113, 147, 97};
    // MKFIVIIFNIICIFP
    int[] query = {131, 128, 147, 113, 99, 113, 113+  147+ 114, 113, 113, 103+57, 113, 147, 97};
    
    //int[] query = {113, 228, 114, 113, 1903};
    ArrayList<Integer> temp = new ArrayList<Integer>();
    for (int mass : query) temp.add(mass);
    
    //int[] query1 = {71, 913, 113, 57, 382, 344, 566};
    ///ArrayList<Integer> temp1 = new ArrayList<Integer>();
    //for (int mass : query1) temp1.add(mass);
    
    /*
    String sub = "ILTVTNFTTCPPDNVSW";
    int cumMass = 0;
    for (int charIndex=0; charIndex<sub.length(); charIndex++) {
      cumMass += AminoAcid.getStandardAminoAcid(sub.charAt(charIndex)).getNominalMass();
    }
    System.out.println("***** MASS " + cumMass);
    */
    
    time = System.currentTimeMillis();
    ArrayList<ArrayList<Integer>> queries = new ArrayList<ArrayList<Integer>>();
    queries.add(temp);
    //queries.add(temp1);
    SequenceGenerator.generateRandomCorrectQueries(sequence, queries, maxItems);
    PSKeywordTree pskt = new PSKeywordTree(queries, sequence, 500);
    System.out.println("-- Done building the key word tree in " + (System.currentTimeMillis() - time)/1000.0 + "s");
    
    //kt.mutationMatch(sequence, results);
    ArrayList<MegaPrefixSuffixMatchObject> results = pskt.collectPrefixSuffixMatches();
    time = System.currentTimeMillis();
    //ArrayList<MatchObject> matches = kt.collectPrefixSuffixMatches();
    //ArrayList<SingleMutationMatchObject> matches = kt.collectMatchesWithOneMutation();
    //System.out.println("-- Done collecting " + matches.size() + " additional prefix-suffix matches in " + (System.currentTimeMillis() - time)/1000.0 + "s");
    
    // check that the matches are correct
    time = System.currentTimeMillis();
    for (MatchObject match : results) {
      System.out.println(match.toString());
      //System.out.println(sequence.getCharAt(2) + " " + sequence.getCharAt(25));
      //System.out.println(match.getMatchAsStringWithFlankingAminoAcids());
      if (!match.isCorrect()) {
        System.out.println(match + " is incorrect");
        System.exit(-9);
      }
    }
    System.out.println("-- Done verifying all queries in " + (System.currentTimeMillis() - time)/1000.0 + "s");
  }
}
