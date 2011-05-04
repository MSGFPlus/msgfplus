package suffixtree.actions;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;

import msgap.Parameters;
import msgap.ScoringParameter;
import msgap.ScoringParameterIterator;
import msgap.results.GappedPeptideResults;
import msgap.results.MSGDResultFileParser;
import sequences.ProteinFastaSequence;
import sequences.ProteinFastaSequences;
import suffixtree.Constants;
import suffixtree.matches.ExactMatchObject;
import suffixtree.trees.KeywordTreeCompact;



public class ExactMatching {

  private static String userHome = System.getProperty("user.home");
  
  
  /**
   * Helper method.
   * @param params
   * @param gpr
   * @param pfs
   * @param out
   */
  private static void matchAndScore(
      Parameters params, 
      GappedPeptideResults gpr, 
      ProteinFastaSequence pfs, 
      PrintWriter out
      ) {

    String matchDir = params.getOutFileName() + "ExactMatchesDir";
    
    int matchedQueries = exactMatching(gpr, pfs, matchDir, false);

    // create a new Scoring Parameter iterator per call
    Iterator<ScoringParameter> sIt = new ScoringParameterIterator(params);
    Scoring.score(params, sIt, gpr, matchDir, matchedQueries, out);
  }
  
  
  
  /**
   * Run the exact matching algorithm and output the results.
   * @param params the parameter file 
   */
  public static void run(Parameters params) {
    
    //long time = System.currentTimeMillis();
    
    // calculate how many results to parse at a time based on the total free memory
    int queryCount = (int)(Runtime.getRuntime().maxMemory() * 0.0006);
    
    System.out.println("\n---> Running " + params.getOutFileName());
    
    PrintWriter out = null;
    try {
      out = new PrintWriter(params.getOutFileName()+".txt");
      // print out the header
      out.println("#"+ExactMatchObject.getSummaryHeader());
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-9);
    }
  
    // The GRC path iterator
    Iterator<GappedPeptideResults> it = new MSGDResultFileParser(params.getGRCPath(), queryCount).iterator();
    
    int spectra = 0, queries = 0;
    while (it.hasNext()) {
      GappedPeptideResults gpr = it.next();
      
      spectra += gpr.getSpectrumCount();
      queries += gpr.getSequenceCount();
      
      if (new File(params.getDBPath()).isDirectory()) {
        ProteinFastaSequences db = new ProteinFastaSequences(params.getDBPath(), sequences.Constants.AMINO_ACIDS_18, params.aaSet(), false); 
        Iterator<ProteinFastaSequence> pfsIt = db.getSequenceIterator();
        boolean append = false;
        String matchDir = params.getOutFileName() + "ExactMatchesDir";
        int matchedQueries = 0;
        while (pfsIt.hasNext()) {
          matchedQueries += exactMatching(gpr, pfsIt.next(), matchDir, append);
          append = true;
        }
        // create a new Scoring Parameter iterator per call
        Iterator<ScoringParameter> sIt = new ScoringParameterIterator(params);
        Scoring.score(params, sIt, gpr, matchDir, matchedQueries, out);
      }
      else {
        ProteinFastaSequence db = new ProteinFastaSequence(params.getDBPath(), sequences.Constants.AMINO_ACIDS_18, params.aaSet());
        matchAndScore(params, gpr, db, out);
      }
      
      System.gc();
    }
    
    out.close(); 
  }
  
  
  
  /**
   * 
   * @param gpr
   * @param db
   * @param matchDir
   * @param append
   * @return
   */
  public static int exactMatching(GappedPeptideResults gpr, 
                                  ProteinFastaSequence db,
                                  String matchDir,
                                  boolean append) {

    System.out.println("\n***** EXACT MATCHING *****");
    long time = System.currentTimeMillis();
    KeywordTreeCompact kt = new KeywordTreeCompact(gpr);
    System.out.printf("Tree built in %d seconds\n", (int)((System.currentTimeMillis()-time)/1000.0));
    int matchedQueries = kt.match(db, matchDir, append);
    kt.closeMatchFiles();
    kt = null;
    return matchedQueries;
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
                    String.format("InputFile %s/Data/Spectra/Asp/ORG014_LTQ_FT_0", userHome)};
    run(new Parameters(in1));
    
    /*
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
    */
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
                    String.format("DBFile %s/Data/Databases/Hsapiens/pro/ipi.HUMAN.v3.78.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Hsapiens/Heck/mzXML/lys", userHome)};
    run(new Parameters(in0));

    String[] in1 = {String.format("OutputFile %s/Data/Spectra/Hsapiens/Heck/mzXML/tryp/results6", userHome), 
                    String.format("DBFile %s/Data/Databases/Hsapiens/pro/ipi.HUMAN.v3.78.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Hsapiens/Heck/mzXML/tryp", userHome)};
    run(new Parameters(in1));
  }
  
  
  public static void runMaize() {
    // Maize
    String[] in0 = {
      String.format("OutputFile %s/Data/Spectra/Maize/SQS02/endosperm-10/output6", userHome), 
      String.format("DBFile %s/Data/Databases/Maize/FASTA", userHome),
      String.format("InputFile %s/Data/Spectra/Maize/SQS02/endosperm-10/", userHome)
      };
    run(new Parameters(in0));
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
    //benchmark();
    runMaize();
  }
  
  
}
