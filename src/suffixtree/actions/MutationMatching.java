package suffixtree.actions;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
//import java.util.ArrayList;
//import java.util.Collections;
//import java.util.HashMap;
//import java.util.HashSet;
import java.util.Iterator;

import msgap.Parameters;
import msgap.ScoringParameter;
import msgap.ScoringParameterIterator;
import msgap.results.GappedPeptideResults;
import msgap.results.MSGDResultFileParser;
//import msutil.AminoAcid;
import sequences.MassSequence;
import sequences.ProteinFastaSequence;
import sequences.ProteinFastaSequences;
import suffixtree.Constants;
//import suffixtree.matches.MatchObject;
import suffixtree.matches.MutMatchObject;
import suffixtree.trees.FRKeywordTreeCompact;



public class MutationMatching {
  
  private static final String userHome = System.getProperty("user.home");

  
  /**
   * Main helper method to do the mutation matching
   * @param params
   * @param gpr
   * @param pfs
   * @param out
   */
  private static void matchAndScore (
      Parameters params, 
      GappedPeptideResults gpr, 
      ProteinFastaSequence pfs, 
      PrintWriter out
      ) {

    String matchDir = params.getOutFileName() + "MutMatchesDir";
    
    int matchedQueries = oneMutationMatching(gpr, pfs, matchDir);

    // create a new Scoring Parameter iterator per call
    Iterator<ScoringParameter> sIt = new ScoringParameterIterator(params);
    Scoring.scoreMutatedMatches(params, pfs, sIt, gpr, matchDir, matchedQueries, out);
  }
  
 
  
  public static void run(Parameters params) {
    
    long time = System.currentTimeMillis();
    
    // calculate how many results to parse at a time based on the total free memory
    int queryCount = (int)(Runtime.getRuntime().maxMemory() * 0.0006);
    
    System.out.println("\n---> Running " + params.getOutFileName());
    
    // open the statistic file
    PrintWriter out = null;
    try {
      out = new PrintWriter(params.getOutFileName()+".txt");
      out.println("#"+MutMatchObject.getSummaryHeader());
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-9);
    }

    // where to get the query files from
    String grcFile = params.getOutFileName() + ".grc";     // the standard suffix added to the mutation run
    
    // the GRC file iterator
    Iterator<GappedPeptideResults> it = new MSGDResultFileParser(grcFile, queryCount).iterator();
    
    int spectra = 0, queries = 0;
    while (it.hasNext()) {
      GappedPeptideResults gpr = it.next();
      
      spectra += gpr.getSpectrumCount();
      queries += gpr.getSequenceCount();
      
      if (new File(params.getDBPath()).isDirectory()) {
        ProteinFastaSequences db =  new ProteinFastaSequences(params.getDBPath(), sequences.Constants.AMINO_ACIDS_18_X, Constants.AA, false);
        Iterator<ProteinFastaSequence> pfsIt = db.getSequenceIterator();
        while (pfsIt.hasNext()) {
          matchAndScore(params, gpr, pfsIt.next(), out);
        }
      }
      else {
        ProteinFastaSequence db = new ProteinFastaSequence(params.getDBPath(), sequences.Constants.AMINO_ACIDS_18_X, Constants.AA);
        matchAndScore(params, gpr, db, out);
      }
     
      out.flush();
      gpr.clear(); // clear the memory
      gpr = null;
      System.gc(); // force garbage collect
    }
    
    time = (System.currentTimeMillis() - time) / 1000;
    System.out.printf("\nElapsed time: %d seconds", time);
    out.close();
  }
  
  
  
  /**
   * One mutation matching and storing the items into a file
   * @param gpr
   * @param db
   * @param matchDir
   * @return
   */
  public static int oneMutationMatching(GappedPeptideResults gpr, 
                                        MassSequence db, 
                                        String matchDir) {

    System.out.println("\n***** MUTATION MATCHING *****");
    FRKeywordTreeCompact kt = new FRKeywordTreeCompact(gpr);

    int matchedQueries = kt.mutationMatch(db, matchDir);
    kt.closeMatchFiles();
    kt = null;
    return matchedQueries;
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
    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_FT_0/results6MutatedR", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translatedR.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Asp/ORG014_LTQ_FT_0", userHome)};
    run(new Parameters(in1));
    
    /*
    String[] in2 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_0/results6Mutated", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Asp/ORG014_LTQ_0", userHome)};
    run(new Parameters(in2));
    
    String[] in3 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_1/results6Mutated", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Asp/ORG014_LTQ_1", userHome)};
    run(new Parameters(in3));

    String[] in4 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_2/results6Mutated", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Asp/ORG014_LTQ_2", userHome)};
    run(new Parameters(in4));
    
    String[] in5 = {String.format("OutputFile %s/Data/Spectra/Asp/ORG014_LTQ_3/results6Mutated", userHome), 
                    String.format("DBFile %s/Data/Databases/Asp/pro/translated.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Asp/ORG014_LTQ_3", userHome)};
    run(new Parameters(in5));
    */
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
  
  public static void bench() {
    // Shewanella
    String[] in0 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT0/benchMutated", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT0", userHome)};
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
    //oli();
    //runAvar();
    //runScerv();
    //runScervMito();
    //runSone();
    bench();
  }
 
}
