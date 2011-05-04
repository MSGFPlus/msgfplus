package suffixtree.actions;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeMap;

import msgap.Parameters;
import msgap.ScoringParameter;
import msgap.ScoringParameterIterator;
import msgap.results.GappedPeptideResults;
import msgap.results.SpectrumMatches;
import msutil.Composition;
import msutil.Peptide;
import sequences.MassSequence;
import suffixtree.matches.ExactMatchObject;
import suffixtree.matches.MatchObject;
import suffixtree.matches.ModMatchObject;
import suffixtree.matches.MutMatchObject;
import suffixtree.misc.ProgressMeter;
import suffixtree.trees.FRKeywordTreeCompact;
import suffixtree.trees.KeywordTreeCompact;



public class Scoring {
  
  private static int readExactMatches(
      GappedPeptideResults gpr, 
      String matchFile, 
      TreeMap<Integer,ArrayList<MatchObject>> matches
      ) {
    
    String matchObjectText;
    int resultsCount = 0;
    try {
      BufferedReader br = new BufferedReader(new FileReader(matchFile));
      while ((matchObjectText = br.readLine()) != null) { // while loop begins here
        ExactMatchObject mo = new ExactMatchObject(gpr.getSequences(), matchObjectText);
        resultsCount++;
        int specId = gpr.getSpecId(mo.getQueryIndex());
        if (!matches.containsKey(specId)) matches.put(specId, new ArrayList<MatchObject>());
        matches.get(specId).add(mo);
      } // end while 
    } // end try
    catch (IOException e) {
      System.err.println("Error: " + e);
    }
    return resultsCount;
  }
  
  
  /**
   * Scores and prints out the modified matches. This implementation prints out
   * the best match(es) with a valid probability
   * @param params the parameter object
   * @param sIt the score parameter iterator
   * @param gpr the gap result object
   * @param out the file to write the results
   * @param totalMatches total number of matches stored in the match directory
   * @param matches the specId to matches map
   */
  public static void score(Parameters params,
                           Iterator<ScoringParameter> sIt,
                           GappedPeptideResults gpr,
                           String matchDir,
                           int totalMatches,
                           PrintWriter out) {

    ProgressMeter pm = new ProgressMeter("\nScoring spectrum-matches", totalMatches, System.out);
    
    int matchFileIndex = 0;
    TreeMap<Integer,ArrayList<MatchObject>> matches = new TreeMap<Integer,ArrayList<MatchObject>>();
    File matchFile = new File(matchDir, String.format(KeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
    readExactMatches(gpr, matchFile.getAbsolutePath(), matches);
    
    //int specCount = 0;
    long time = System.currentTimeMillis();
    int scoredMatches = 0;

    while (sIt.hasNext()) {
      ScoringParameter sp = sIt.next();
      int specId = sp.getSpecID();
      
      if (specId > matches.lastKey()) {
        matches.clear();
        matchFile = new File(matchDir, String.format(KeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
        if (!matchFile.exists()) break;
        //System.out.println(matchFile.getAbsolutePath());
        readExactMatches(gpr, matchFile.getAbsolutePath(), matches);
      }
        
      //boolean hasMatch = false;
      if (matches.containsKey(specId)) {
        // score all the matches
        SpectrumMatches sm = new SpectrumMatches(sp, params);
        
        float bestProb = sm.addMatches(matches.get(specId));
        scoredMatches += matches.get(specId).size();
        pm.update(scoredMatches);
        
        if (bestProb < 1.0f) {
          
          for (MatchObject mo : sm.getMatches()) {
            
            if (mo.getProb()<=bestProb) {
              // offset is the ParentMass - Theoretical Mass
              float offset = sp.getOriginalParentMass()-mo.getPeptide().getMass()-(float)Composition.H2O;
              String line = mo.getSummaryLine(gpr.getFileName(specId), gpr.getScanNumber(specId), gpr.getActmethod(specId), gpr.getPrecursorMass(specId), gpr.getCharge(specId), offset);
              out.println(line);
            }
          }
        }
       
      }
      if (specId >= gpr.getMaxId()) {
        if (specId > gpr.getMaxId()) {
          System.out.println("WARNING: Potential spectrum that was not scored");
        }
        break; // we are done iterating the iterator
      }
    }
    System.out.println();
    
    time = System.currentTimeMillis() - time;
    System.out.printf("Scored %d matches in %d seconds.\n", scoredMatches, time/1000);
  }
  
  
  

  private static int readModMatches(MassSequence db, 
                                    GappedPeptideResults gpr, 
                                    String matchFile, 
                                    TreeMap<Integer,ArrayList<MatchObject>> matches) {
    String matchObjectText;
    int resultsCount = 0;
    try {
      BufferedReader br = new BufferedReader(new FileReader(matchFile));
      while ((matchObjectText = br.readLine()) != null) { // while loop begins here
        ModMatchObject mo = new ModMatchObject(db, gpr.getSequences(), matchObjectText);
        resultsCount++;
        int specId = gpr.getSpecId(mo.getQueryIndex());
        if (!matches.containsKey(specId)) matches.put(specId, new ArrayList<MatchObject>());
        matches.get(specId).add(mo);
      } // end while 
    } // end try
    catch (IOException e) {
      System.err.println("Error: " + e);
    }
    return resultsCount;
  }
  
  
  private static int readMutMatches (
      MassSequence db, 
      GappedPeptideResults gpr, 
      String matchFile, 
      TreeMap<Integer,ArrayList<MatchObject>> matches
      ) {

    String matchObjectText;
    int resultsCount = 0;
    try {
      BufferedReader br = new BufferedReader(new FileReader(matchFile));
      while ((matchObjectText = br.readLine()) != null) { // while loop begins here
        MutMatchObject mo = new MutMatchObject(db, gpr.getSequences(), matchObjectText);
        resultsCount++;
        int specId = gpr.getSpecId(mo.getQueryIndex());
        if (!matches.containsKey(specId)) matches.put(specId, new ArrayList<MatchObject>());
        matches.get(specId).add(mo);
      } // end while 
    } // end try
    catch (IOException e) {
      System.err.println("Error: " + e);
    }
    return resultsCount;
  }
  
  
  
  /**
   * Scores and prints out the mutated matches. This implementation prints out
   * the best match(es) with a valid probability
   * @param params the parameter object
   * @param sIt the score parameter iterator
   * @param gpr the gap result object
   * @param out the file to write the results
   * @param totalMatches total number of matches stored in the match directory
   * @param matches the specId to matches map
   */
  public static void scoreMutatedMatches(Parameters params,
                                         MassSequence db,
                                         Iterator<ScoringParameter> sIt,
                                         GappedPeptideResults gpr,
                                         String matchDir,
                                         int totalMatches,
                                         PrintWriter out) {

    ProgressMeter pm = new ProgressMeter("\nScoring spectrum-matches", totalMatches, System.out);
    
    int matchFileIndex = 0;
    TreeMap<Integer,ArrayList<MatchObject>> matches = new TreeMap<Integer,ArrayList<MatchObject>>();
    File matchFile = new File(matchDir, String.format(FRKeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
    readMutMatches(db, gpr, matchFile.getAbsolutePath(), matches);
    
    //int specCount = 0;
    long time = System.currentTimeMillis();
    int scoredMatches = 0;

    while (sIt.hasNext()) {
      ScoringParameter sp = sIt.next();
      int specId = sp.getSpecID();
      
      if (specId > matches.lastKey()) {
        matches.clear();
        matchFile = new File(matchDir, String.format(FRKeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
        //System.out.println(matchFile.getAbsolutePath());
        readMutMatches(db, gpr, matchFile.getAbsolutePath(), matches);
      }
        
      //boolean hasMatch = false;
      if (matches.containsKey(specId)) {
        // score all the matches
        SpectrumMatches sm = new SpectrumMatches(sp, params);
        
        float bestProb = sm.addMatches(matches.get(specId));
        scoredMatches += matches.get(specId).size();
        pm.update(scoredMatches);
        
        if (bestProb < 1.0f) {
          
          for (MatchObject mo : sm.getMatches()) {
            
            if (mo.getProb()<=bestProb) {
              // offset is the ParentMass - Theoretical Mass
              float offset = sp.getOriginalParentMass()-mo.getPeptide().getMass()-(float)Composition.H2O;
              String line = mo.getSummaryLine(gpr.getFileName(specId), gpr.getScanNumber(specId), gpr.getActmethod(specId), gpr.getPrecursorMass(specId), gpr.getCharge(specId), offset);
              out.println(line);
            }
          }
        }
       
      }
      if (specId >= gpr.getMaxId()) {
        if (specId > gpr.getMaxId()) {
          System.out.println("WARNING: Potential spectrum that was not scored");
        }
        break; // we are done iterating the iterator
      }
    }
    System.out.println();
    
    time = System.currentTimeMillis() - time;
    System.out.printf("Scored %d matches in %d seconds.\n", scoredMatches, time/1000);
  }
  
  
  
  /**
   * Scores and prints out the modified matches. This implementation prints out
   * the best match(es) with a valid probability
   * @param params the parameter object
   * @param sIt the score parameter iterator
   * @param gpr the gap result object
   * @param out the file to write the results
   * @param totalMatches total number of matches stored in the match directory
   * @param matches the specId to matches map
   */
  public static void scoreModdedMatches(Parameters params,
                                        MassSequence db,
                                        Iterator<ScoringParameter> sIt,
                                        GappedPeptideResults gpr,
                                        String matchDir,
                                        int totalMatches,
                                        PrintWriter out) {

    ProgressMeter pm = new ProgressMeter("\nScoring spectrum-matches", totalMatches, System.out);
    
    int matchFileIndex = 0;
    TreeMap<Integer,ArrayList<MatchObject>> matches = new TreeMap<Integer,ArrayList<MatchObject>>();
    File matchFile = new File(matchDir, String.format(FRKeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
    readModMatches(db, gpr, matchFile.getAbsolutePath(), matches);
    
    //int specCount = 0;
    long time = System.currentTimeMillis();
    int scoredMatches = 0;

    while (sIt.hasNext()) {
      ScoringParameter sp = sIt.next();
      int specId = sp.getSpecID();
      
      if (specId > matches.lastKey()) {
        matches.clear();
        matchFile = new File(matchDir, String.format(FRKeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
        //System.out.println(matchFile.getAbsolutePath());
        readModMatches(db, gpr, matchFile.getAbsolutePath(), matches);
      }
        
      //boolean hasMatch = false;
      if (matches.containsKey(specId)) {
        // score all the matches
        SpectrumMatches sm = new SpectrumMatches(sp, params);
        
        float bestProb = sm.addMatches(matches.get(specId));
        scoredMatches += matches.get(specId).size();
        pm.update(scoredMatches);
        
        if (bestProb < 1.0f) {
          
          for (MatchObject mo : sm.getMatches()) {
            
            if (mo.getProb()<=bestProb) {
              // offset is the ParentMass - Theoretical Mass
              float offset = sp.getOriginalParentMass()-mo.getPeptide().getMass()-(float)Composition.H2O;
              String line = mo.getSummaryLine(gpr.getFileName(specId), gpr.getScanNumber(specId), gpr.getActmethod(specId), gpr.getPrecursorMass(specId), gpr.getCharge(specId), offset);
              out.println(line);
            }
          }
        }
       
      }
      if (specId >= gpr.getMaxId()) {
        if (specId > gpr.getMaxId()) {
          System.out.println("WARNING: Potential spectrum that was not scored");
        }
        break; // we are done iterating the iterator
      }
    }
    System.out.println();
    
    time = System.currentTimeMillis() - time;
    System.out.printf("Scored %d matches in %d seconds.\n", scoredMatches, time/1000);
  }
  
  

  /**
   * Scores and prints out the matches. This implementation prints out
   * the best match(es) with a valid probability
   * @param params the parameter object
   * @param sIt the score parameter iterator
   * @param gpr the gap result object
   * @param out the file to write the results
   * @param matches the specId to matches map
   */
  public static void score(Parameters params,
                           Iterator<ScoringParameter> sIt,
                           GappedPeptideResults gpr,
                           PrintWriter out, 
                           HashMap<Integer,ArrayList<MatchObject>> matches) {

    ProgressMeter pm = new ProgressMeter("\nScoring spectrum-matches", matches.size(), System.out);
    
    int specCount = 0;
    long time = System.currentTimeMillis();
    int scoredMatches = 0;

    while (sIt.hasNext()) {
      ScoringParameter sp = sIt.next();
      int specId = sp.getSpecID();
      
      //boolean hasMatch = false;
      if (matches.containsKey(specId)) {
        // score all the matches
        SpectrumMatches sm = new SpectrumMatches(sp, params);
        pm.update(++specCount);
        
        float bestProb = sm.addMatches(matches.get(specId));
        scoredMatches += matches.get(specId).size();
        
        if (bestProb < 1.0f) {
          
          for (MatchObject mo : sm.getMatches()) {
            
            if (mo.getProb()<=bestProb) {
              // offset is the ParentMass - Theoretical Mass
              float offset = sp.getOriginalParentMass()-mo.getPeptide().getMass()-(float)Composition.H2O;
              String line = mo.getSummaryLine(gpr.getFileName(specId), gpr.getScanNumber(specId), gpr.getActmethod(specId), gpr.getPrecursorMass(specId), gpr.getCharge(specId), offset);
              out.println(line);
            }
          }
        }
       
      }
      if (specId >= gpr.getMaxId()) {
        if (specId > gpr.getMaxId()) {
          System.out.println("WARNING: Potential spectrum that was not scored");
        }
        break; // we are done iterating the iterator
      }
    }
    System.out.println();
    
    time = System.currentTimeMillis() - time;
    System.out.printf("Scored %d matches in %d seconds.\n", scoredMatches, time/1000);
  }
  
  
  
  // score a particular spectrum
  public static void score(Parameters params, String filename, int scanNum) {
    Iterator<ScoringParameter> i = new ScoringParameterIterator(params);
    
    String peptideF = "R.GLGILDSALNELQGDTLDGETVFK.L";
    String peptide = "GLGILDSALNELQGDTLDGETVFK";
    while (i.hasNext()) {
      
      ScoringParameter sp = i.next();
      String[] tokens = sp.getSpecFileName().split("/");
      String name = tokens[tokens.length-1].trim();
      //System.out.println(sp.getScanNum());
      //System.out.println(name);
      if (name.equals(filename) && scanNum==sp.getScanNum()) {
        System.out.println(name); 
        SpectrumMatches sm = new SpectrumMatches(sp, params);
        Peptide pep = new Peptide(peptide);
        int score = sm.getScore(pep);
        float prob = sm.getProbability(pep, peptideF, score);
        System.out.printf("%s\t%d\t%.3e\n", name, score, prob);
      }
      
    }
  }
  
  
  public static void main(String[] args) {
    String userHome = System.getProperty("user.home");
    
    String[] in3 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT3/output6", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT3", userHome)};
    score(new Parameters(in3), "ShewFed037_LTQFT_1_12Nov04_Pegasus_0804-4_dta.ms2", 10998);  
  }
  
  
  
}
