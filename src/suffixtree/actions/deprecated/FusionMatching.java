package suffixtree.actions.deprecated;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import msgap.Parameters;
//import msgap.ScoringParameter;
//import msgap.ScoringParameterIterator;
import msgap.results.GappedPeptideResults;
import msgap.results.MSGDResultFileParser;
import msutil.AminoAcidSet;
import sequences.ProteinFastaSequence;
import sequences.ProteinFastaSequences;
import suffixtree.Constants;
//import suffixtree.actions.Scoring;
import suffixtree.matches.MatchObject;
import suffixtree.matches.MegaPrefixSuffixMatchObject;
import suffixtree.misc.MatchingStats;
import suffixtree.trees.PSKeywordTree;

public class FusionMatching {

  
  private static HashSet<Integer> matchAndScore(Parameters params, 
                                                GappedPeptideResults gpr, 
                                                ProteinFastaSequence pfs, 
                                                PrintWriter out,
                                                PrintWriter stats,
                                                MatchingStats statObj) {
    
    float probCutOff = params.getSpecProb();
    
    // create the result
    HashMap<Integer,ArrayList<MatchObject>> matches = new HashMap<Integer,ArrayList<MatchObject>>();
    
  	fusionMatching(gpr, pfs, matches, stats);
  	
  	// create a new Scoring Parameter iterator per call
  	//Iterator<ScoringParameter> sIt = new ScoringParameterIterator(params);
  	//Scoring.score(sIt, params, gpr, out, probCutOff, statObj, matches);
  	
    // create the spectra that did not have a satisfactory match score
    HashSet<Integer> matchedSpecs = new HashSet<Integer>();
    for (int specId : matches.keySet()) {
      boolean passed = false;
      for (MatchObject mo : matches.get(specId)) {
        if (mo.getProb()<=probCutOff) {
           passed = true;
           break;
        }
      }
      if (passed) matchedSpecs.add(specId);
    }
    
    return matchedSpecs;
  }
    
  
  
  /**
   * Run the exact matching algorithm and output the results.
   * @param params the parameter file 
   */
  public static void run(Parameters params) {
    
    // open the statistic file
    PrintWriter stats = null;
    PrintWriter out = null;
    PrintWriter grc = null;
    try {
      out = new PrintWriter(params.getOutFileName()+"PrefixSuffix.txt");
      stats = new PrintWriter(params.getOutFileName()+"PrefixSuffix.log");
      grc = new PrintWriter(params.getOutFileName()+"PrefixSuffix.grc");
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-9);
    }
    
    // statistics object
    MatchingStats statObj = new MatchingStats();
   
    String grcFile = params.getOutFileName() + "Splicing.grc";     // the standard suffix added to the ps run
    Iterator<GappedPeptideResults> it = new MSGDResultFileParser(grcFile, 1000000).iterator();
    while (it.hasNext()) {
      GappedPeptideResults gpr = it.next();
 
      HashSet<Integer> matchedSpecs = new HashSet<Integer>();
      
      // fusion matching
      if (new File(params.getDBPath()).isDirectory()) {
        ProteinFastaSequences db = new ProteinFastaSequences(params.getDBPath(), sequences.Constants.AMINO_ACIDS_18, params.aaSet(), false); 
        Iterator<ProteinFastaSequence> pfsIt = db.getSequenceIterator();
        while (pfsIt.hasNext()) {
          matchedSpecs.addAll(matchAndScore(params, gpr, pfsIt.next(), out, stats, statObj));
        }
      }
      else {
        ProteinFastaSequence db = new ProteinFastaSequence(params.getDBPath(), sequences.Constants.AMINO_ACIDS_18, params.aaSet());
        matchedSpecs.addAll(matchAndScore(params, gpr, db, out, stats, statObj));
      }
   
      GappedPeptideResults leftovers = gpr.generateGPR(matchedSpecs);
      leftovers.toFile(grc);
      out.flush();
      gpr.clear();
    }
    stats.close(); out.close(); grc.close();
  }
  
  
  
  /**
   * 
   * @param gpr
   * @param db
   * @param matches
   * @param stats
   * @return
   */
  public static void fusionMatching(GappedPeptideResults gpr, 
                                    ProteinFastaSequence db, 
                                    HashMap<Integer,ArrayList<MatchObject>> matches,
                                    PrintWriter stats) {

    System.out.println("\n***** FUSION MATCHING *****");
    long time = System.currentTimeMillis();
    PSKeywordTree pskt = new PSKeywordTree(gpr.getSequences(), db, 10);
    
    if (stats!=null) {
      stats.println("----- FUSION MATCHING -----");
      stats.println("build keyword tree:" + (System.currentTimeMillis()-time)/1000.0);
    }
    
    //time = System.currentTimeMillis();
    ArrayList<MegaPrefixSuffixMatchObject> results = pskt.collectPrefixSuffixMatches();
    if (stats!=null) stats.println("searching keyword tree:" + (System.currentTimeMillis()-time)/1000.0);
    
    // populate the matches data-structure
    HashSet<Integer> seenIds = new HashSet<Integer>();
    HashSet<Integer> matchedQueries = new HashSet<Integer>();
    for (MatchObject mo : results) {
      int specId = gpr.getSpecId(mo.getQueryIndex());
      if (!matches.containsKey(specId)) matches.put(specId, new ArrayList<MatchObject>());
      matches.get(specId).add(mo);
      matchedQueries.add(mo.getQueryIndex());
      seenIds.add(specId);
    }
    
    System.out.println("\n-- Summary --");
    System.out.printf("Matched %d out of %d spectra\n", matches.size(), gpr.getSpectrumCount());
    System.out.printf("Matched %d out of %d queries.\n", results.size(), gpr.getSequenceCount());
    System.out.printf("Average of %.2f matches per spectrum\n", results.size()/(float)matches.size());
    System.out.printf("Average of %.2f matches per query\n", results.size()/(float)matchedQueries.size());
    System.out.printf("Done querying in %.2f seconds\n", (System.currentTimeMillis()-time)/1000.0);
    
    if (stats!=null) {
      stats.println("total spectra:" + gpr.getSpectrumCount());
      stats.println("total identified spectra:" + matches.size());
      stats.println("total queries (gapped sequences):" + gpr.getSequenceCount());
      stats.println("total matched queries (gapped sequences):" + results.size());
    }
  }
  
  
  public static void batchRun() {
    String userHome = System.getProperty("user.home");
    
    
    // Asp
    /*
    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_FT_0/allOutput", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in1));
    
    String[] in2 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_0/allOutput", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in2));
    
    String[] in3 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_1/allOutput", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in3));

    String[] in4 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_2/allOutput", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in4));
    
    String[] in5 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_3/allOutput", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in5));

    
    // Csp
    String[] in6 = {String.format("OutputFile %s/Data/Spectra/Csp/ORG033_LTQ_Orb_0/allOutput", userHome), 
                    String.format("DBFile %s/Data/Databases/Csp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in6));

    String[] in7 = {String.format("OutputFile %s/Data/Spectra/Csp/ORG033_LTQ_0/allOutput", userHome), 
                    String.format("DBFile %s/Data/Databases/Csp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in7));

    
    // Avar
    String[] in8 = {String.format("OutputFile %s/Data/Spectra/Avar/ORG013_LTQ_0/allOutput", userHome), 
                    String.format("DBFile %s/Data/Databases/Avar/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in8));

    String[] in9 = {String.format("OutputFile %s/Data/Spectra/Avar/ORG013_LTQ_1/allOutput", userHome), 
                    String.format("DBFile %s/Data/Databases/Avar/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
   run(new Parameters(in9));
   */
    
   // Scerv
   String[] in10 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LTQ_Orb_0/allOutput", userHome), 
                    String.format("DBFile %s/Data/Databases/Scerv/pro", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
   run(new Parameters(in10));

   String[] in11 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LTQ_FT_0/allOutput", userHome), 
                    String.format("DBFile %s/Data/Databases/Scerv/pro", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
   run(new Parameters(in11));
   
   String[] in12 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LTQ_0/allOutput", userHome), 
                    String.format("DBFile %s/Data/Databases/Scerv/pro", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
   run(new Parameters(in12));

   String[] in13 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LCQ_0/allOutput", userHome), 
                    String.format("DBFile %s/Data/Databases/Scerv/pro", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
   run(new Parameters(in13));
  }
  
  
  public static void run() {
    String userHome = System.getProperty("user.home");
    
    // Ecoli
    /*
    String[] in = {String.format("OutputFile %s/Data/Spectra/Ecoli/ORG045_LTQ_Orb_0/allOutput", userHome), 
                   String.format("DBFile %s/Data/Databases/Ecoli/pro/translated.fasta", userHome),
                   String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    */
    
    // Yeast
    String[] in = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LTQ_Orb_0/results", userHome), 
                   String.format("DBFile %s/Data/Databases/Scerv/pro", userHome),
                   String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    
    Parameters params = new Parameters(in);
    run(params);
  }
  
  
  public static void runHuman() {
    String userHome = System.getProperty("user.home");
    
    // Human
    String[] in = {String.format("OutputFile %s/Data/Spectra/Hsapiens/Heck/mzXML/lys/results", userHome), 
                   String.format("DBFile %s/Data/Databases/Hsapiens/translated", userHome),
                   String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    
    Parameters params = new Parameters(in);
    run(params);
  }
  
  
  public static void main(String[] args) {
    Constants.AA = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
    //batchRun();
    //run();
    runHuman();
  }
  
  
}
