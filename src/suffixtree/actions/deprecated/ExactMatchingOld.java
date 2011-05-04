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
import sequences.ProteinFastaSequence;
import sequences.ProteinFastaSequences;
import suffixtree.Constants;
//import suffixtree.actions.Scoring;
import suffixtree.matches.ExactMatchObject;
import suffixtree.matches.MatchObject;
import suffixtree.misc.MatchingStats;
import suffixtree.trees.KeywordTree;



public class ExactMatchingOld {

  private static String userHome = System.getProperty("user.home");
  
  private static HashSet<Integer> matchAndScore(Parameters params, 
                                                GappedPeptideResults gpr, 
                                                ProteinFastaSequence pfs, 
                                                PrintWriter out,
                                                PrintWriter stats,
                                                MatchingStats statObj) {

    float probCutOff = params.getSpecProb();

    // create the result
    HashMap<Integer,ArrayList<MatchObject>> matches = new HashMap<Integer,ArrayList<MatchObject>>();

    exactMatching(gpr, pfs, matches);

    // create a new Scoring Parameter iterator per call
    //Iterator<ScoringParameter> sIt = new ScoringParameterIterator(params);
    System.out.println();
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
    
    long time = System.currentTimeMillis();
    
    // calculate how many results to parse at a time based on the total free memory
    int queryCount = (int)(Runtime.getRuntime().maxMemory() * 0.0006);
    
    System.out.println("\n---> Running " + params.getOutFileName());
    
    // open the statistic file
    PrintWriter stats = null;
    PrintWriter out = null;
    PrintWriter grc = null;
    try {
      out = new PrintWriter(params.getOutFileName()+".txt");
      stats = new PrintWriter(params.getOutFileName()+"Exact.log");
      grc = new PrintWriter(params.getOutFileName()+"Exact.grc");
      
      // print out the header
      out.print("# SpecFilePath\t");
      out.print("ScanNum\t");
      out.print("ActivationMethod\t");
      out.print("PrecursorMz\t");
      out.print("Charge\t");
      out.print(ExactMatchObject.getSummaryHeader()+"\t");
      out.println("Delta");
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-9);
    }
  
    // statistics object
    MatchingStats statObj = new MatchingStats();

    // The GRC path iterator
    Iterator<GappedPeptideResults> it = new MSGDResultFileParser(params.getGRCPath(), queryCount).iterator();
    
    int spectra = 0, queries = 0;
    while (it.hasNext()) {
      GappedPeptideResults gpr = it.next();
      
      spectra += gpr.getSpectrumCount();
      queries += gpr.getSequenceCount();
      
      HashSet<Integer> matchedSpecs = new HashSet<Integer>();
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
      
      System.gc();
    }
    
    // print out the statistics
    statObj.printMatchStats(stats, 1e-10f);
    statObj.printMatchStats(stats, 1e-11f);
    statObj.printMatchStats(stats, 1e-12f);
    statObj.printMatchStats(stats, 1e-13f);
    statObj.printMatchStats(stats, 1e-14f);
    statObj.printOutScoreDist(stats);
    stats.printf("%d spectra with %d queries in the grp files.\n", spectra, queries);
    time = (System.currentTimeMillis() - time) / 1000;
    stats.printf("\nElapsed time: %d seconds\n", time);
    
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
  public static void exactMatching(GappedPeptideResults gpr, 
                                   ProteinFastaSequence db, 
                                   HashMap<Integer,ArrayList<MatchObject>> matches) {

    System.out.println("\n***** EXACT MATCHING *****");
    long time = System.currentTimeMillis();
    KeywordTree kt = new KeywordTree(gpr.getSequences());
    System.out.printf("Tree built in %d seconds\n", (int)((System.currentTimeMillis()-time)/1000.0));
    
    ArrayList<ExactMatchObject> results = new ArrayList<ExactMatchObject>();
    time = System.currentTimeMillis();
    kt.match(db, results);
    
    // populate the matches data-structure
    HashSet<Integer> seenIds = new HashSet<Integer>();
    HashSet<Integer> matchedQueries = new HashSet<Integer>();
    for (ExactMatchObject mo : results) {
      int specId = gpr.getSpecId(mo.getQueryIndex());
      if (!matches.containsKey(specId)) matches.put(specId, new ArrayList<MatchObject>());
      matches.get(specId).add(mo);
      matchedQueries.add(mo.getQueryIndex());
      seenIds.add(specId);
    }
    System.out.printf("Done querying in %d seconds\n", (int)((System.currentTimeMillis()-time)/1000.0));
    
    System.out.println("\n-- Summary --");
    System.out.printf("Matched %d out of %d spectra\n", matches.size(), gpr.getSpectrumCount());
    System.out.printf("Matched %d out of %d queries.\n", results.size(), gpr.getSequenceCount());
    System.out.printf("Average of %.2f matches per spectrum\n", results.size()/(float)matches.size());
    System.out.printf("Average of %.2f matches per query\n", results.size()/(float)matchedQueries.size());
    System.out.printf("Done querying in %.2f seconds\n", (System.currentTimeMillis()-time)/1000.0);
  }
  
  
  public static void runAvar() {
    // Avar
    String[] in0 = {String.format("OutputFile %s/Data/Spectra/Avar/ORG013_LTQ_0/results6", userHome), 
                    String.format("DBFile %s/Data/Databases/Avar/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in0));

    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Avar/ORG013_LTQ_1/results6", userHome), 
                    String.format("DBFile %s/Data/Databases/Avar/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
   run(new Parameters(in1));
  }
  
  
  public static void runCsp() {
    // Csp
    String[] in0 = {String.format("OutputFile %s/Data/Spectra/Csp/ORG033_LTQ_Orb_0/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Csp/pro/translated.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Csp/ORG033_LTQ_Orb_0", userHome), 
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in0));

    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Csp/ORG033_LTQ_0/results", userHome), 
                    String.format("DBFile %s/Data/Databases/Csp/pro/translated.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Csp/ORG033_LTQ_0", userHome), 
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in1));
  }
  
  public static void runSone() {
    // Shewanella
    String[] in0 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT0/output6", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    //String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants_R.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT0", userHome)};
    run(new Parameters(in0));

    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT1/output6", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    //String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants_R.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT1", userHome)};
    run(new Parameters(in1));
    
    String[] in2 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT2/output6", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    //String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants_R.fasta", userHome),            
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT2", userHome)};
    run(new Parameters(in2));

    String[] in3 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT3/output6", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    //String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants_R.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT3", userHome)};
    run(new Parameters(in3));

    String[] in4 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT4/output6", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    //String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants_R.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT4", userHome)};
    run(new Parameters(in4));

    String[] in5 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT5/output6", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    //String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants_R.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT5", userHome)};
    run(new Parameters(in5));

  }
  
  
  
  public static void benchmark() {
    String[] in = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT0/bench", userHome), 
                   //String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants_R.fasta", userHome),
                   String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                   String.format("InputFile %s/Data/Spectra/Sone/LTQFT0", userHome)};
    run(new Parameters(in));
  }
  
  
  
  public static void runAsp() {
    // Asp
    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_FT_0/results6", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in1));
    
    String[] in2 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_0/results6", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in2));
    
    String[] in3 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_1/results6", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in3));

    String[] in4 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_2/results6", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in4));
    
    String[] in5 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_3/results6", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in5));
  }
  
  
  
  public static void runEcoli() {
    // Ecoli 
    String[] in0 = {String.format("OutputFile %s/Data/Spectra/Ecoli/ORG045_LTQ_Orb_0/output", userHome), 
                    String.format("DBFile %s/Data/Databases/Ecoli/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in0));

    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Ecoli/ORG045_LTQ_Orb_1/output", userHome), 
                    String.format("DBFile %s/Data/Databases/Ecoli/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in1));
    
    String[] in2 = {String.format("OutputFile %s/Data/Spectra/Ecoli/ORG045_LTQ_Orb_2/output", userHome), 
                    String.format("DBFile %s/Data/Databases/Ecoli/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in2));

    String[] in3 = {String.format("OutputFile %s/Data/Spectra/Ecoli/ORG045_LTQ_Orb_3/output", userHome), 
                    String.format("DBFile %s/Data/Databases/Ecoli/pro/translated.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in3));
  }
  
  
  public static void runScerv() {
    // Scerv
    String[] in0 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LTQ_Orb_0/results6", userHome), 
                    //String.format("DBFile %s/Data/Databases/Scerv/pro", userHome),
                    String.format("DBFile %s/Data/Databases/Scerv/orf/orf_trans_all.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in0));

    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LTQ_FT_0/results6", userHome), 
                    //String.format("DBFile %s/Data/Databases/Scerv/pro", userHome),
                    String.format("DBFile %s/Data/Databases/Scerv/orf/orf_trans_all.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in1));
    
    String[] in2 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LTQ_0/results6", userHome), 
                    //String.format("DBFile %s/Data/Databases/Scerv/pro", userHome),
                    String.format("DBFile %s/Data/Databases/Scerv/orf/orf_trans_all.fasta", userHome), 
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in2));

    String[] in3 = {String.format("OutputFile %s/Data/Spectra/Scerv/ORG105_LCQ_0/results6", userHome), 
                    //String.format("DBFile %s/Data/Databases/Scerv/pro", userHome),
                    String.format("DBFile %s/Data/Databases/Scerv/orf/orf_trans_all.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in3));
  }
  
  
  public static void runHsapiens() {
    // Human
    String[] in0 = {String.format("OutputFile %s/Data/Spectra/Hsapiens/Heck/mzXML/lys/results6", userHome), 
                    String.format("DBFile %s/Data/Databases/Hsapiens/pro/proteins.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in0));

    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Hsapiens/Heck/mzXML/tryp/results6", userHome), 
                    String.format("DBFile %s/Data/Databases/Hsapiens/pro/proteins.fasta", userHome),
                    String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    run(new Parameters(in1));
  }
  
  
  
  public static void run() {
    String userHome = System.getProperty("user.home");
    String[] in = {String.format("OutputFile %s/Data/Spectra/Csp/ORG033_LTQ_Orb_0/results", userHome), 
                   String.format("DBFile %s/Data/Databases/Csp/pro/translated.fasta", userHome),
                   String.format("SpecProb %e", Constants.PROB_CUTOFF)};
    
    Parameters params = new Parameters(in);
    run(params);
  }
  
  
  public static void main(String[] args) {
    //run();
    //runScerv();
    //runAsp();
    //runCsp();
    //runEcoli();
    //runAvar();
    //runHsapiens();
    //runSone();
    benchmark();
  }
  
  
}
