package suffixtree.actions.deprecated;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;

import sequences.MassSequence;
import sequences.ProteinFastaSequence;
import suffixtree.matches.MatchObject;
import msgap.results.GappedPeptideResults;
import msgap.results.MSGDResultFileParser;


/**
 * Given a set of spectra run the exact matching and mutation matching on the
 * remaining items
 * @author jung
 *
 */
public class FullAnalysis {

  
  /*
  public static void run3steps(String generalName, String specFile, String queryFile, String dbFile, String outFile, String statsFile) {
    
    float mutProbCutOff = 1e-13f;
    
    // open the protein sequence file
    ProteinFastaSequence db = new ProteinFastaSequence(dbFile);
    
    // open the statistic file
    PrintWriter stats = null;
    try {
      stats = new PrintWriter(statsFile);
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-9);
    }
    
    // parse the query file
    GappedPeptideResults gpr = new MSGDResultFileParser(queryFile).iterator().next();
    
    // create the result
    HashMap<Integer,ArrayList<MatchObject>> matches = new HashMap<Integer,ArrayList<MatchObject>>();
    
    // matching exactly
    ExactMatching.exactMatching(gpr, db, matches, stats);
    
    // Score the matches
    Scoring.score(specFile, outFile, stats, suffixtree.Constants.PROB_CUTOFF, matches);
    
    // create the spectra that did not have a satisfactory match score
    HashSet<Integer> matchedSpecs = new HashSet<Integer>();
    for (int specId : matches.keySet()) {
      boolean passed = false;
      for (MatchObject mo : matches.get(specId)) {
        if (mo.getProb()<=suffixtree.Constants.PROB_CUTOFF) {
           passed = true;
           break;
        }
      }
      if (passed) matchedSpecs.add(specId);
    }
    
    //GappedPeptideResults leftovers = gpr.removeSpecs(matchedSpecs);
    //System.out.println("Unmmatched spectra " + leftovers.getSpectrumCount());

    
    
    TreeSet<Integer> additionalIds = new TreeSet<Integer>();
    
    // run it 2 more times for the correct parent masses
    matches.clear();
    String inputFile = System.getProperty("user.home")+String.format("/Data/Testground/mutations/%splus.grc", generalName);
    gpr = new MSGDResultFileParser(inputFile).iterator().next();
    ExactMatching.exactMatching(gpr, db, matches, stats);
    
    // Score the matches
    Scoring.score(System.getProperty("user.home")+"/Data/Testground/mutations/merged10000plus.ms2", outFile, stats, suffixtree.Constants.PROB_CUTOFF, matches);
    
    // create the spectra that did not have a satisfactory match score
    for (int specId : matches.keySet()) {
      boolean passed = false;
      for (MatchObject mo : matches.get(specId)) {
        if (mo.getProb()<=suffixtree.Constants.PROB_CUTOFF) {
           passed = true;
           break;
        }
      }
      if (passed) {
        matchedSpecs.add(specId);
        additionalIds.add(specId);
      }
    }
    
    //leftovers = gpr.removeSpecs(matchedSpecs);
    //System.out.println("Unmmatched spectra " + leftovers.getSpectrumCount());

    
    
    matches.clear();
    inputFile = System.getProperty("user.home")+String.format("/Data/Testground/mutations/%sminus.grc", generalName);
    gpr = new MSGDResultFileParser(inputFile).iterator().next();
    ExactMatching.exactMatching(gpr, db, matches, stats);
    
    // Score the matches
    Scoring.score(System.getProperty("user.home")+"/Data/Testground/mutations/merged10000minus.ms2", outFile, stats, suffixtree.Constants.PROB_CUTOFF, matches);
    
    // create the spectra that did not have a satisfactory match score
    for (int specId : matches.keySet()) {
      boolean passed = false;
      for (MatchObject mo : matches.get(specId)) {
        if (mo.getProb()<=suffixtree.Constants.PROB_CUTOFF) {
           passed = true;
           break;
        }
      }
      if (passed) {
        matchedSpecs.add(specId);
        additionalIds.add(specId);
      }
    }

    // parse the query file, re-parse the query, ffs
    gpr = new MSGDResultFileParser(queryFile).iterator().next();
    
    GappedPeptideResults leftovers = gpr.removeSpecs(matchedSpecs);
    System.out.println("Unmmatched spectra " + leftovers.getSpectrumCount());
    System.out.println("Matched spectra " + matchedSpecs.size());

    try {
      PrintWriter pw = new PrintWriter("/home/jung/Desktop/additionalIds.txt");
      for (int specId : additionalIds) {
        pw.println(specId);
      }
      pw.close();
    }
    catch (IOException ioe) {
      System.err.println(ioe);
    }    
    
    
    
    
    matches.clear();
    MutationMatching.oneMutationMatching(leftovers, db, matches, stats);
    Scoring.score(specFile, outFile, stats, mutProbCutOff, matches);
    
    //MutationMatching.printMutationTable(matches, mutProbCutOff, db);
    MutationMatching.collectRepeatingPointMutations(matches, mutProbCutOff, db);
  }
  */
  
  private static void printProteins(MassSequence db, String outfile, float probCutOff, HashMap<Integer,ArrayList<MatchObject>> matches) {
    HashMap<String,ArrayList<String>> groupedMatches = new HashMap<String,ArrayList<String>>();

    // group the matches by proteins, which is the protein name
    for (int specId : matches.keySet()) {
      for (MatchObject mo : matches.get(specId)) {
        if (mo.getProb()<=probCutOff) {
          String annotation = db.getAnnotation(mo.getStart());
          if (!groupedMatches.containsKey(annotation)) {
            groupedMatches.put(annotation, new ArrayList<String>());
          }
          groupedMatches.get(annotation).add(mo.toString());
        }
      }
    }
    
    // We use a key which puts the count first to facilitate sorting
    TreeSet<String> sortedIds = new TreeSet<String>();
    for (String key : groupedMatches.keySet()) {
      sortedIds.add(String.format("%05d", groupedMatches.get(key).size()) + "_" + key);
    }
    
    // sort everything in reverse order by the match count
    try {
      PrintWriter pw = new PrintWriter(outfile);
      for (String codedString : sortedIds.descendingSet()) {
        String[] tokens = codedString.split("_", 2);
        pw.println(tokens[1] + "\t" + groupedMatches.get(tokens[1]).size());
        for (String match : groupedMatches.get(tokens[1])) {
          pw.println(match);
        }
        pw.println();
      }
      pw.close();
    }
    catch (IOException ioe) {
      System.err.println(ioe);
    }
  }
  
  
  public static void runWithRescoring(String queryFile, ProteinFastaSequence db, String outFilePrefix) {
    
    float mutProbCutOff = 1e-13f;
    
    // open the statistic and output file
    PrintWriter stats = null;    // open the output file
    PrintWriter out = null;
    PrintWriter mutatedOut = null;
    try {
      out = new PrintWriter(outFilePrefix+"Exact.txt");
      stats = new PrintWriter(outFilePrefix + ".log");
      mutatedOut = new PrintWriter(outFilePrefix+"Mutated.txt");
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-9);
    }
    
    // parse the query file
    GappedPeptideResults gpr = new MSGDResultFileParser(queryFile).iterator().next();
    
    // create the result
    HashMap<Integer,ArrayList<MatchObject>> matches = new HashMap<Integer,ArrayList<MatchObject>>();
    
    // matching exactly
    ExactMatchingOld.exactMatching(gpr, db, matches);
    
    // Score the matches
    //Scoring.score(gpr, out, stats, suffixtree.Constants.PROB_CUTOFF, matches);
    
    // create the spectra that did not have a satisfactory match score
    HashSet<Integer> matchedSpecs = new HashSet<Integer>();
    for (int specId : matches.keySet()) {
      boolean passed = false;
      for (MatchObject mo : matches.get(specId)) {
        if (mo.getProb()<=suffixtree.Constants.PROB_CUTOFF) {
           passed = true;
           break;
        }
      }
      if (passed) matchedSpecs.add(specId);
    }
    GappedPeptideResults leftovers = gpr.generateGPR(matchedSpecs);
    System.out.println("Unmmatched spectra " + leftovers.getSpectrumCount());
    
    matches.clear();
    //MutationMatching.oneMutationMatching(leftovers, db, matches, stats);
    //Scoring.score(specFile, outFilePrefix+"Mutated.txt", stats, leftovers, mutProbCutOff, matches);
    //Scoring.scoreMutation(leftovers, mutatedOut, stats, mutProbCutOff, matches);
    
    // print out the mutated results as proteins
    printProteins(db, outFilePrefix+"Proteins.txt", mutProbCutOff, matches);
    
    //MutationMatching.printMutationTable(matches, outFilePrefix+"MutationTable.txt", mutProbCutOff, db);
    //MutationMatching.collectRepeatingPointMutations(matches, mutProbCutOff, db);
    
    stats.close(); out.close(); mutatedOut.close();
  }
  
  
  
  public static void run(String queryFile, String dbFile, String outFilePrefix) {
    
    //float mutProbCutOff = 1e-13f;
    
    // open the protein sequence file
    ProteinFastaSequence db = new ProteinFastaSequence(dbFile);
    
    // open the statistic file
    PrintWriter stats = null;
    PrintWriter exactOut = null;
    PrintWriter mutatedOut = null;
    try {
      stats = new PrintWriter(outFilePrefix + ".log");
      exactOut = new PrintWriter(outFilePrefix + "Exact.txt");
      mutatedOut = new PrintWriter(outFilePrefix+"Mutated.txt");
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-9);
    }
    
    // parse the query file
    GappedPeptideResults gpr = new MSGDResultFileParser(queryFile).iterator().next();
    
    // create the result
    HashMap<Integer,ArrayList<MatchObject>> matches = new HashMap<Integer,ArrayList<MatchObject>>();
    
    // matching exactly
    ExactMatchingOld.exactMatching(gpr, db, matches);
    
    // Score the matches
    //Scoring.score(gpr, exactOut, stats, suffixtree.Constants.PROB_CUTOFF, matches);
    
    // create the spectra that did not have a satisfactory match score
    HashSet<Integer> matchedSpecs = new HashSet<Integer>();
    for (int specId : matches.keySet()) {
      boolean passed = false;
      for (MatchObject mo : matches.get(specId)) {
        if (mo.getProb()<=suffixtree.Constants.PROB_CUTOFF) {
           passed = true;
           break;
        }
      }
      if (passed) matchedSpecs.add(specId);
    }
    GappedPeptideResults leftovers = gpr.generateGPR(matchedSpecs);
    System.out.println("Unmmatched spectra " + leftovers.getSpectrumCount());
    
    matches.clear();
    //MutationMatching.oneMutationMatching(leftovers, db, matches, stats);
    //Scoring.score(leftovers, mutatedOut, stats, mutProbCutOff, matches);
    
    //MutationMatching.printMutationTable(matches, outFilePrefix+"MutationTable.txt", mutProbCutOff, db);
    //MutationMatching.collectRepeatingPointMutations(matches, mutProbCutOff, db);
    
    stats.close();
    mutatedOut.close();
    exactOut.close();
  }


  
  public static void runMouse() {
    String generalName = "MmusHeartMito1000";
    
    String userHome = System.getProperty("user.home");
    String queryFile = userHome + String.format("/Data/Gap/spectra/%s.grc", generalName);
    String outFilePrefix = userHome + String.format("/Data/Gap/spectra/%s", generalName);
    
    // databases
    //String dbFile = userHome + "/Data/Databases/Mmus/pro/mtp.fasta";
    String dbFile = userHome + "/Data/Databases/Mmus/pro/Mus_musculus.NCBIM37.59.pep.all.fasta";
    //String dbFile = userHome + "/Data/Gap/spectra/MmusHeartMitoIds.fasta";
    
    //String rDbFile = userHome + "/Data/Databases/ShewDB/SOne_uniprot_plus_contaminants_R.fasta";
    ProteinFastaSequence db = new ProteinFastaSequence(dbFile);
    //MassSequence rDb = new ProteinFastaSequence(rDbFile);
   
    //run(specFile, queryFile, dbFile, outFile, statsFile);
    //runWithRescoring(specFile, queryFile, rDb, outFilePrefix);
    runWithRescoring(queryFile, db, outFilePrefix);
  }
  
  
  public static void runRat() {
    String generalName = "MmusHeartMito";
    
    String userHome = System.getProperty("user.home");
    String queryFile = userHome + String.format("/Data/Gap/spectra/%s.grc", generalName);
    String outFilePrefix = userHome + String.format("/Data/Gap/spectra/%sRat", generalName);
    
    // databases
    //String dbFile = userHome + "/Data/Databases/Mmus/pro/mtp.fasta";
    //String dbFile = userHome + "/Data/Databases/Mmus/pro/Mus_musculus.NCBIM37.59.pep.all.fasta";
    String dbFile = userHome + "/Data/Databases/Rnor/pro/Rattus_norvegicus.RGSC3.4.59.pep.all.fasta";
    //String dbFile = userHome + "/Data/Gap/spectra/MmusHeartMitoIds.fasta";
    
    //String rDbFile = userHome + "/Data/Databases/ShewDB/SOne_uniprot_plus_contaminants_R.fasta";
    ProteinFastaSequence db = new ProteinFastaSequence(dbFile);
    //MassSequence rDb = new ProteinFastaSequence(rDbFile);
   
    //run(specFile, queryFile, dbFile, outFile, statsFile);
    //runWithRescoring(specFile, queryFile, rDb, outFilePrefix);
    runWithRescoring(queryFile, db, outFilePrefix);
  }
  
  
  
  
  public static void main(String[] args) {
    //runYeast();
    //runShew();
    //runRat();
    //runMouse();
  }
  
  
}
