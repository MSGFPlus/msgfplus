package suffixtree.actions;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import msgap.Parameters;
import msgap.ScoringParameter;
import msgap.ScoringParameterIterator;
import msgap.results.GappedPeptideResults;
import msgap.results.MSGDResultFileParser;
import sequences.MassSequence;
import sequences.ProteinFastaSequence;
import sequences.ProteinFastaSequences;
import suffixtree.Constants;
import suffixtree.matches.MatchObject;
import suffixtree.matches.ModMatchObject;
import suffixtree.trees.FRKeywordTree;
import suffixtree.trees.FRKeywordTreeCompact;



public class ModificationMatching {
  
  private static final String userHome = System.getProperty("user.home");
  
  private static void matchAndScore(Parameters params,
                                    GappedPeptideResults gpr, 
                                    ProteinFastaSequence pfs, 
                                    PrintWriter out) {

    String matchDir = params.getOutFileName() + "ModMatchesDir";
    int matchedQueries = blindMatching(gpr, pfs, matchDir);
    
    // create a new Scoring Parameter iterator per call
    Iterator<ScoringParameter> sIt = new ScoringParameterIterator(params);
    Scoring.scoreModdedMatches(params, pfs, sIt, gpr, matchDir, matchedQueries, out);
  }
  
  
  public static void run(Parameters params) {
    
    long time = System.currentTimeMillis();
    
    // calculate how many results to parse at a time based on the total free memory
    int queryCount = (int)(Runtime.getRuntime().maxMemory() * 0.0005);
    
    System.out.println("\n---> Running " + params.getOutFileName());
    
    PrintWriter out = null;
    try {
      out = new PrintWriter(params.getOutFileName()+".txt");
      out.println("#"+ModMatchObject.getSummaryHeader());
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-9);
    }
    
    String grcFile = params.getOutFileName() + ".grc";     // the standard suffix added to the mutation run
    
    // the GRC file iterator
    Iterator<GappedPeptideResults> it = new MSGDResultFileParser(grcFile, queryCount, 6).iterator();
    
    int spectra = 0, queries = 0;
    while (it.hasNext()) {
      GappedPeptideResults gpr = it.next();
      
      spectra += gpr.getSpectrumCount();
      queries += gpr.getSequenceCount();
      
      if (new File(params.getDBPath()).isDirectory()) {
        ProteinFastaSequences db = new ProteinFastaSequences(params.getDBPath(), sequences.Constants.AMINO_ACIDS_18, params.aaSet(), false); 
        Iterator<ProteinFastaSequence> pfsIt = db.getSequenceIterator();
        while (pfsIt.hasNext()) {
          matchAndScore(params, gpr, pfsIt.next(), out);
        }
      }
      else {
        ProteinFastaSequence db = new ProteinFastaSequence(params.getDBPath(), sequences.Constants.AMINO_ACIDS_18, params.aaSet());
        matchAndScore(params, gpr, db, out);
      }
      
      // generate the 
      out.flush();
      
      gpr.clear(); // clear the memory
      gpr = null;
      
      System.gc();
    }
    
    time = (System.currentTimeMillis() - time) / 1000;
    System.out.printf("\nElapsed time: %d seconds", time);
    out.close();
  }

  
  
  /**
   * Alternative matching that stores the matches in a file temporarily 
   * @param gpr the gapped peptide queries
   * @param db the database to match
   * @param matchDir the directory to store the matches
   * @return the number of matched queries
   */
  public static int blindMatching(GappedPeptideResults gpr, 
                                  MassSequence db, 
                                  String matchDir) {

    System.out.println("\n***** BLIND MATCHING *****");
    //long time = System.currentTimeMillis();
    FRKeywordTreeCompact kt = new FRKeywordTreeCompact(gpr);
    
    int matchedQueries = kt.blindMatch(db, matchDir);
    kt.closeMatchFiles();     // close the match file and read the results
    kt = null;
    return matchedQueries;
    
    /*
    String matchObjectText;
    int resultsCount = 0;
    HashSet<Integer> seenIds = new HashSet<Integer>();
    HashSet<Integer> matchedQueries = new HashSet<Integer>();
    try {
      BufferedReader br = new BufferedReader(new FileReader(matchFile));
      while ((matchObjectText = br.readLine()) != null) { // while loop begins here
        ModMatchObject mo = new ModMatchObject(db, gpr.getSequences(), matchObjectText);
        resultsCount++;
        int specId = gpr.getSpecId(mo.getQueryIndex());
        if (!matches.containsKey(specId)) matches.put(specId, new ArrayList<MatchObject>());
        matches.get(specId).add(mo);
        matchedQueries.add(mo.getQueryIndex());
        seenIds.add(specId);
      } // end while 
    } // end try
    catch (IOException e) {
      System.err.println("Error: " + e);
    }
    
    
    System.out.println("\n-- Summary --");
    System.out.printf("Matched %d out of %d spectra\n", matches.size(), gpr.getSpectrumCount());
    System.out.printf("Matches %d out of %d queries.\n", resultsCount, gpr.getSequenceCount());
    System.out.printf("Average of %.2f matches per spectrum\n", resultsCount/(float)matches.size());
    System.out.printf("Average of %.2f matches per query\n", resultsCount/(float)matchedQueries.size());
    System.out.printf("Done querying in %.2f seconds\n", (System.currentTimeMillis()-time)/1000.0);
    
    return gpr.generateGPR(seenIds);
    */
  }
  
  
  public static GappedPeptideResults blindMatching(GappedPeptideResults gpr, 
                                                   MassSequence db, 
                                                   HashMap<Integer,ArrayList<MatchObject>> matches,
                                                   PrintWriter stats) {

    System.out.println("\n***** BLIND MATCHING *****");
    long time = System.currentTimeMillis();
    FRKeywordTree kt = new FRKeywordTree(gpr.getSequences());

    if (stats!=null) {
      stats.println("----- BLIND MATCHING -----");
      stats.println("build keyword tree:" + (System.currentTimeMillis()-time)/1000.0);
    }
    
    ArrayList<ModMatchObject> results = new ArrayList<ModMatchObject>();
    time = System.currentTimeMillis();
    kt.blindMatch(db, results);
    if (stats!=null) stats.println("searching keyword tree:" + (System.currentTimeMillis()-time)/1000.0);

    // populate the matches data-structure
    HashSet<Integer> seenIds = new HashSet<Integer>();
    HashSet<Integer> matchedQueries = new HashSet<Integer>();
    for (ModMatchObject mo : results) {
      int specId = gpr.getSpecId(mo.getQueryIndex());
      if (!matches.containsKey(specId)) matches.put(specId, new ArrayList<MatchObject>());
      matches.get(specId).add(mo);
      matchedQueries.add(mo.getQueryIndex());
      seenIds.add(specId);
    }

    System.out.println("\n-- Summary --");
    System.out.printf("Matched %d out of %d spectra\n", matches.size(), gpr.getSpectrumCount());
    System.out.printf("Matches %d out of %d queries.\n", results.size(), gpr.getSequenceCount());
    System.out.printf("Average of %.2f matches per spectrum\n", results.size()/(float)matches.size());
    System.out.printf("Average of %.2f matches per query\n", results.size()/(float)matchedQueries.size());
    System.out.printf("Done querying in %.2f seconds\n", (System.currentTimeMillis()-time)/1000.0);
    
    // print the matches into a file
    if (stats!=null) {
      stats.println("total spectra:" + gpr.getSpectrumCount());
      stats.println("total identified spectra:" + matches.size());
      stats.println("total queries (gapped sequences):" + gpr.getSequenceCount());
      stats.println("total matched queries (gapped sequences):" + results.size());
    }
    return gpr.generateGPR(seenIds);

  }
  
  
  
  public static void runAvar() {
    // Avar
    String[] in0 = {String.format("OutputFile %s/Data/Spectra/Avar/ORG013_LTQ_0/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Avar/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in0));

    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Avar/ORG013_LTQ_1/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Avar/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
   run(new Parameters(in1));
  }
  
  public static void runCsp() {
    // Csp
    String[] in0 = {String.format("OutputFile %s/Data/Spectra/Csp/ORG033_LTQ_Orb_0/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Csp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in0));

    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Csp/ORG033_LTQ_0/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Csp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in1));
  }
  
  
  public static void runAsp() {
    // Asp
    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_FT_0/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in1));
    
    String[] in2 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_0/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in2));
    
    String[] in3 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_1/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in3));

    String[] in4 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_2/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in4));
    
    String[] in5 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_3/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in5));
  }
  
  
  public static void runEcoli() {
    // Ecoli 
    String[] in3 = {String.format("OutputFile %s/Data/Spectra/Ecoli/ORG045_LTQ_Orb_3/output", userHome), 
                    String.format("DBFile %s/Data/Databases/Ecoli/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in3));

    String[] in2 = {String.format("OutputFile %s/Data/Spectra/Ecoli/ORG045_LTQ_Orb_2/output", userHome), 
                    String.format("DBFile %s/Data/Databases/Ecoli/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in2));

    String[] in0 = {String.format("OutputFile %s/Data/Spectra/Ecoli/ORG045_LTQ_Orb_0/output", userHome), 
                    String.format("DBFile %s/Data/Databases/Ecoli/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in0));

    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Ecoli/ORG045_LTQ_Orb_1/output", userHome), 
                    String.format("DBFile %s/Data/Databases/Ecoli/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in1));
  }
  
  
  public static void bench() {
    // Shewanella
    String[] in0 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT0/benchModded", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT0", userHome)};
    run(new Parameters(in0));
  }
  

  public static void runHsapiens() {
    // Human IPI run
    String[] in0 = {String.format("OutputFile %s/Data/Spectra/Hsapiens/Heck/mzXML/tryp/output6Modded", userHome), 
                    String.format("DBFile %s/Data/Databases/Hsapiens/pro/ipi.HUMAN.v3.78.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Hsapiens/Heck/mzXML/tryp", userHome)};
    run(new Parameters(in0));
  }
  
  
  public static void runSone() {
    // Shewanella
    String[] in0 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT0/output6", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT0", userHome)};
    run(new Parameters(in0));

    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT1/output6", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT1", userHome)};
    run(new Parameters(in1));
    
    String[] in2 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT2/output6", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT2", userHome)};
    run(new Parameters(in2));

    String[] in3 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT3/output6", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT3", userHome)};
    run(new Parameters(in3));

    String[] in4 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT4/output6", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT4", userHome)};
    run(new Parameters(in4));

    String[] in5 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT5/output6", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT5", userHome)};
    run(new Parameters(in5));
  }
    
  
  public static void runScerv() {
    // Scerv
    String[] in0 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LTQ_Orb_0/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Scerv/pro", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in0));

    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LTQ_FT_0/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Scerv/pro", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in1));
    
    String[] in2 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LTQ_0/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Scerv/pro", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in2));

    String[] in3 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LCQ_0/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Scerv/pro", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in3));
  }
  
  public static void runScervMito() {
    // Scerv
    String[] in0 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LTQ_Orb_0/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Scerv/pro/chrmtp.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in0));

    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LTQ_FT_0/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Scerv/pro/chrmtp.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in1));
    
    String[] in2 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LTQ_0/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Scerv/pro/chrmtp.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in2));

    String[] in3 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LCQ_0/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Scerv/pro/chrmtp.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in3));
  }
  
  public static void run() {
    String userHome = System.getProperty("user.home");

    /*
    String[] in = {String.format("OutputFile %s/Data/Spectra/Ecoli/ORG045_LTQ_7/allResultsExact", userHome), 
                   String.format("DBFile %s/Data/Databases/Ecoli/pro/translated.fasta", userHome),
                   String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    */
    
    String[] in = {String.format("OutputFile %s/Data/Spectra/Csp/ORG033_LTQ_Orb_0/results", userHome), 
                   String.format("DBFile %s/Data/Databases/Csp/pro/translated.fasta", userHome),
                   String.format("SpecProb %e", Constants.PROB_CUTOFF)};

    Parameters params = new Parameters(in);
    run(params);
  }
  
  
  public static void main(String[] args) {
    //run();
    //runAsp();
    //runCsp();
    //runEcoli();
    //runAvar();
    //runScerv();
    //runScervMito();
    //runSone();
    bench();
    //runHsapiens();
  }
 
}
