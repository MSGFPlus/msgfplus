package suffixtree.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

import sequences.ProteinFastaSequence;
import suffixtree.matches.ExactMatchObject;
import suffixtree.trees.KeywordTree;

/**
 * Main class to match a set of gap peptides to the a given data base and output
 * the matches
 * @author jung
 *
 */
public class GapMatcher {

  /**
   * Main function to search the database given a result file from msgap.
   * @param inFile the path to the result file in msgap
   * @param database the path to the fasta database
   * @return the map object that associates a list of MatchObject's to the spectrum scan number
   */
  public static HashMap<Integer,ArrayList<ExactMatchObject>> matchDB(String inFile, String database) {

    long time = System.currentTimeMillis();
    HashMap<Integer,String> ids = new HashMap<Integer,String>();
    HashMap<Integer,ArrayList<Integer>> ids2sequences = new HashMap<Integer,ArrayList<Integer>>();
    HashMap<Integer,Integer> sequences2ids = new HashMap<Integer,Integer>();
    ArrayList<ArrayList<Integer>> sequences = new ArrayList<ArrayList<Integer>>();
    BufferedReader fid;
    try {
      fid = new BufferedReader(new FileReader(inFile));
      String line, ident = null;
      int currentSpec = -1;
      while((line=fid.readLine())!=null) {
        if (line.charAt(0)=='#') {
          // parse the identification line
          StringTokenizer tk = new StringTokenizer(line);
          int index = Integer.parseInt(tk.nextToken().substring(1));
          ident = tk.nextToken();
          currentSpec = index;
          ids.put(index, ident);
          ids2sequences.put(index, new ArrayList<Integer>());
        }
        else {
          // parse the array
          String trimmed = line.trim();
          StringTokenizer tk = new StringTokenizer(trimmed.substring(1, trimmed.length()-1), ", ");
          ArrayList<Integer> sequence = new ArrayList<Integer>();
          while (tk.hasMoreTokens()) {
            int mass = Integer.parseInt(tk.nextToken());
            //if (mass>500) System.out.println("Fail mass " + mass);
            sequence.add(mass);
          }
          ids2sequences.get(currentSpec).add(sequences.size());
          sequences2ids.put(sequences.size(), currentSpec);
          sequences.add(sequence);
        }
      }
      fid.close();
    }
    catch(IOException ioe) {
      System.err.println(ioe.getStackTrace());
      System.exit(-9);
    }
    
    System.out.println("Sequences " + sequences2ids.size());
    System.out.println("Specs " + ids.size());

    // Create the SuffixTree
    KeywordTree kt = new KeywordTree(sequences);
    ProteinFastaSequence db = new ProteinFastaSequence(database);
    ArrayList<ExactMatchObject> results = new ArrayList<ExactMatchObject>();
    kt.match(db, results);
    
    // group the matches by spec id
    HashMap<Integer,ArrayList<ExactMatchObject>> matches = new HashMap<Integer,ArrayList<ExactMatchObject>>();
    for (ExactMatchObject mo : results) {
      int specId = sequences2ids.get(mo.getQueryIndex());
      if (!matches.containsKey(specId)) matches.put(specId, new ArrayList<ExactMatchObject>());
      matches.get(specId).add(mo);
    }

    System.out.println("Elapsed time: " + (System.currentTimeMillis()-time)/1000 + " second");
    return matches;
  }
  
  /**
   * This function can be called to print the results of the match into a file
   * @param outFile the output file
   * @param matches the matches from the matching function
   */
  public static void writeToFile(String outFile, HashMap<Integer,ArrayList<ExactMatchObject>> matches) {
    BufferedWriter fout;
    try {
      fout = new BufferedWriter(new FileWriter(outFile));
      for (int specId : matches.keySet()) {
        fout.write(String.format("#%d\n", specId));
        for (ExactMatchObject match : matches.get(specId)) {
          fout.write(match.getMatchAsString()+"\n");
        }
      }
      fout.close();
    }
    catch (IOException ioe) {
      System.err.println(ioe.getStackTrace());
      System.exit(-9);
    }
  }
  
  public static void main(String[] args) {
    String userHome = System.getProperty("user.home");
    String results = userHome+"/Desktop/msgap.txt";
    String db = userHome+"/Data/Databases/ShewDB/wholeGenomeShewTranslationWithContam.fasta";
    String outFile = userHome+"/Desktop/matches.txt";
    
    writeToFile(outFile, matchDB(results, db));
  }
}
