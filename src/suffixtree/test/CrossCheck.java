package suffixtree.test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;

import msutil.AminoAcid;

import sequences.ProteinFastaSequence;
import suffixtree.Constants;
import suffixtree.matches.ExactMatchObject;
import suffixtree.trees.HashedIntegerGappedSuffixTree;
import suffixtree.trees.KeywordTree;

/**
 * This class is to check that both Keyword Tree and Suffix Tree produce the
 * same results
 * @author jung
 *
 */
public class CrossCheck {

  
  public static void generateRandomQueries(ArrayList<ArrayList<Integer>> queries, int count) {
    Random r = new Random();
    
    AminoAcid[] stdAa = AminoAcid.getStandardAminoAcids(); 
    
    for (int i=0; i<count; i++) {
      ArrayList<Integer> query = new ArrayList<Integer>();
      int baseMass = Constants.MIN_QUERY_MASS+Constants.MAX_GAP_MASS+1;
      int maxMass = baseMass+r.nextInt(Constants.MAX_QUERY_MASS-baseMass);
      int totalMass = 0;
      while(true) {
        int maxGapMass = r.nextInt(Constants.MAX_GAP_MASS+1);
        int cumMass = stdAa[r.nextInt(stdAa.length)].getNominalMass();
        while (true) {
          int nextMass = stdAa[r.nextInt(stdAa.length)].getNominalMass();
          if (cumMass+nextMass > maxGapMass) {
            break;
          }
          cumMass += nextMass;
        }
        if (totalMass + cumMass > maxMass) {
          break;
        }
        totalMass += cumMass;
        query.add(cumMass);
      }
      queries.add(query);
    }
  }
  
  public static void printQuery(ArrayList<Integer> query) {
    for (int mass : query) {
      System.err.print(mass + ", ");
    }
    System.err.println();
  }
  
  public static void main(String[] args) {
    String userHome = System.getProperty("user.home");
    String fastaFile;
    
    int correctQueries = 100000;
    int incorrectQueries = correctQueries*9;
    
    fastaFile = userHome+"/Data/Databases/test.fasta";
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
    System.out.println("--- Number of characters in fasta file: " + sequence.getSize());
    
    ArrayList<ArrayList<Integer>> queries = new ArrayList<ArrayList<Integer>>();
    SequenceGenerator.generateRandomCorrectQueries(sequence, queries, correctQueries);
    generateRandomQueries(queries, incorrectQueries);

    time = System.currentTimeMillis();
    HashedIntegerGappedSuffixTree higst = new HashedIntegerGappedSuffixTree(sequence);
    System.out.println("--- Done building the tree in " + (System.currentTimeMillis() - time)/1000.0 + "s");
    
    time = System.currentTimeMillis();
    KeywordTree kt = new KeywordTree(queries);
    System.out.println("--- Done building the keyword tree in " + (System.currentTimeMillis() - time)/1000.0 + "s");
    
    System.out.println("-- Total queries " + queries.size());
    
    ArrayList<ExactMatchObject> ktResults = new ArrayList<ExactMatchObject>();
    time = System.currentTimeMillis();
    kt.match(sequence, ktResults);
    System.out.println("--- Keyword feeding in " + (System.currentTimeMillis() - time)/1000.0 + "s");
    HashSet<Long> coors = new HashSet<Long>();
    for (ExactMatchObject mo : ktResults) {
      coors.add(((long)mo.getStart())<<32 | mo.getEnd());
    }

    time = System.currentTimeMillis();
    HashSet<ExactMatchObject> results = new HashSet<ExactMatchObject>();
    for (ArrayList<Integer> query : queries) {
      higst.search(query, results);
    }
    System.out.println("--- Gapped Suffix Tree querying in " + (System.currentTimeMillis() - time)/1000.0 + "s");
    
    // check for correctness
    for (ArrayList<Integer> query : queries) {
      results.clear();
      higst.search(query, results);
      for (ExactMatchObject mo : results) {
        long key = ((long)mo.getStart())<<32 | mo.getEnd();
        //System.out.println(mo);
        if (!coors.contains(key)) {
          printQuery(query);
          System.err.println(mo + " not found");
          System.exit(-9);
        }
      }
    }
    System.out.println("--- Done cross-checking all answers seem to match");
    
  }
  
}
