package suffixtree.test.deprecated;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.TreeMap;

import msgap.NewMSGappedDictionary;
import msgap.Parameters;
import msgap.results.GappedPeptideResults;
import msgap.results.MSGDResultFileParser;
import msgap.results.SpectrumMatches;
import msutil.AminoAcid;
import msutil.SpectraIterator;
import msutil.Spectrum;
import parser.MS2SpectrumParser;

import sequences.ProteinFastaSequence;
import suffixtree.actions.deprecated.ExactMatchingOld;
import suffixtree.matches.ExactMatchObject;
import suffixtree.matches.MatchObject;
import suffixtree.matches.MutMatchObject;
import suffixtree.trees.HashedIntegerGappedSuffixTree;
import suffixtree.trees.KeywordTree;



/**
 * The idea of this class is to analyze a directory of gap peptide results. 
 * For example, exact match, prefix-suffix match, mutation match, etc.
 * @author jung
 *
 */
public class Analysis {

  private static class Param {
    String input;
    String output;
    String db;
    String mutatedOutput;
  }
   
  public static void runVerify(String dirPath, String dbPath, String exactResultFile, String mutatedResultFile) {
    ProteinFastaSequence db = new ProteinFastaSequence(dbPath);
        
    GappedPeptideResults gpr = new MSGDResultFileParser(dirPath).iterator().next();
    
    // do exact matching and discard those items that do have match
    GappedPeptideResults forMutationSearch = matchingExperiment(gpr, db, exactResultFile);
    System.out.println("Sequences with no matches in the DB " + forMutationSearch.getSequenceCount()+"\n");
    
    //oneMutationMatching(forMutationSearch, db, mutatedResultFile);
  }
  
  /**
   * Collects the statistics to verify that the mutation search is working.
   * @param gpr
   * @param db
   * @param outPath
   * @return
   */
  public static GappedPeptideResults matchingExperiment(GappedPeptideResults gpr, ProteinFastaSequence db, String outPath) {
    
    long time = System.currentTimeMillis();
    System.out.println("***** Starting to search for exact matches ");
    KeywordTree kt = new KeywordTree(gpr.getSequences());
    ArrayList<ExactMatchObject> results = new ArrayList<ExactMatchObject>();
    kt.match(db, results);
    
    HashSet<Integer> seenIds = new HashSet<Integer>();
    // maps the specId to the array of matches
    HashMap<Integer,ArrayList<ExactMatchObject>> matches = new HashMap<Integer,ArrayList<ExactMatchObject>>();
    for (ExactMatchObject mo : results) {
      int specId = gpr.getSpecId(mo.getQueryIndex());
      if (!matches.containsKey(specId)) matches.put(specId, new ArrayList<ExactMatchObject>());
      matches.get(specId).add(mo);
      seenIds.add(specId);
    }
    
    // print out the results
    try {
      PrintWriter fout = new PrintWriter(outPath);
      for (int specId : matches.keySet()) {
        for (ExactMatchObject emo : matches.get(specId)) {
          fout.println(emo.toString());
        }
      }
      fout.close();
    }
    catch (IOException e) {}
    
    System.out.printf("Total specs with matches %d out of %d spec queries. Total matches in the DB %d\n", matches.size(), gpr.getSpectrumCount(), gpr.getSequenceCount());
    System.out.printf("On average there are %.2f db matches per spectrum\n", gpr.getSequenceCount()/(float)matches.size());
    System.out.println("Done quering in " + ((System.currentTimeMillis()-time) / 1000.0) + " seconds");
    return gpr.generateGPR(seenIds);
  }
  
  
  public static GappedPeptideResults exactMatching(GappedPeptideResults gpr, 
                                                   HashedIntegerGappedSuffixTree higst, 
                                                   HashMap<Integer,ArrayList<ExactMatchObject>> matches,
                                                   PrintWriter stats) {

    System.out.println("***** Starting to search for exact matches using the suffix tree");
    long time = System.currentTimeMillis();
    ArrayList<ExactMatchObject> results = new ArrayList<ExactMatchObject>();
    HashSet<ExactMatchObject> tempResults = new HashSet<ExactMatchObject>();
    for (int queryIndex = 0; queryIndex < gpr.getSequenceCount(); queryIndex++) {
      higst.search(gpr.getSequenceAt(queryIndex), tempResults);
      for (ExactMatchObject emo : tempResults) {
        emo.setQueryIndex(queryIndex);
        results.add(emo);
      }
      tempResults.clear();
    }
    if (stats!=null) stats.println("searching suffix tree:" + (System.currentTimeMillis()-time)/1000.0);
    
    HashSet<Integer> seenIds = new HashSet<Integer>();
    // maps the specId to the array of matches
    //HashMap<Integer,ArrayList<ExactMatchObject>> matches = new HashMap<Integer,ArrayList<ExactMatchObject>>();
    int totalMatches = 0;
    for (ExactMatchObject mo : results) {
    int specId = gpr.getSpecId(mo.getQueryIndex());
    if (!matches.containsKey(specId)) matches.put(specId, new ArrayList<ExactMatchObject>());
    matches.get(specId).add(mo);
    totalMatches++;
    seenIds.add(specId);
    }
    
    System.out.printf("Total specs with matches %d out of %d total spec. Total matches in the DB %d out of %d total queries\n", matches.size(), gpr.getSpectrumCount(), totalMatches, gpr.getSequenceCount());
    System.out.printf("On average there are %.2f db matches per spectrum\n", totalMatches/(float)matches.size());
    System.out.println("Done querying in " + ((System.currentTimeMillis()-time) / 1000.0) + " seconds");
    
    if (stats!=null) {
    stats.println("total spectra:" + gpr.getSpectrumCount());
    stats.println("total identified spectra:" + matches.size());
    stats.println("total queries (gapped sequences):" + gpr.getSequenceCount());
    stats.println("total matched queries (gapped sequences):" + totalMatches);
    }
    return gpr.generateGPR(seenIds);
  }

  


  

  
  
  public static void runSearch(String dirPath, String dbPath, String exactResultFile, String mutatedResultFile) {
    ProteinFastaSequence db = new ProteinFastaSequence(dbPath);
        
    GappedPeptideResults gpr = new MSGDResultFileParser(dirPath).iterator().next();
    
    HashMap<Integer,ArrayList<MatchObject>> matches = new HashMap<Integer,ArrayList<MatchObject>>();
    
    // do exact matching and discard those items that do have match
    ExactMatchingOld.exactMatching(gpr, db, matches);
    
    // print out the results
    try {
      PrintWriter fout = new PrintWriter(exactResultFile);
      for (int specId : matches.keySet()) {
        for (MatchObject emo : matches.get(specId)) {
          fout.println(emo.toString());
        }
      }
      fout.close();
    }
    catch (IOException e) {
      System.err.println(e);
    }
    
    
    //HashMap<Integer,ArrayList<SingleMutationMatchObject>> mutatedMatches = new HashMap<Integer,ArrayList<SingleMutationMatchObject>>();
    //oneMutationMatching(forMutationSearch, db, mutatedMatches, mutatedResultFile);
  }


  /**
   * This routine benchmarks the speed of the SuffixTree vs the KeywordTree in 
   * matching gapped patterns.
   * @param p The parameter object
   * @param db The path to the database
   * @param statsFile the path to the statistic files.
   */
  public static void compareKeywordSuffixTree(Parameters p, String db, String statsFile) {
    
    ProteinFastaSequence sequences = new ProteinFastaSequence(db);
    PrintWriter stats = null;
    try {
      stats = new PrintWriter(statsFile);
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-9);
    }
    
    stats.printf("Database size %d\n", sequences.getSize());
    // create the suffix tree
    //long time = System.currentTimeMillis();
    //HashedIntegerGappedSuffixTree higst = new HashedIntegerGappedSuffixTree(sequences);
    //stats.println("\n***** suffix tree:" + db);
    //stats.println("build suffix tree:" + (System.currentTimeMillis()-time)/1000.0);
    
    // assume the input is an MS2 file for now
    for (String msFile : p.specFiles()) {
      
      String resultFile = p.getOutFileName() + ".grc";
      System.out.println("Processing " + resultFile);
      
      GappedPeptideResults gpr = new MSGDResultFileParser(resultFile).iterator().next();
      
      // maps the sequence id to the list of matches
      HashMap<Integer,ArrayList<MatchObject>> matches = new HashMap<Integer,ArrayList<MatchObject>>();
      
      //stats.println("\n***** suffix tree tree:" + msFile);
      //exactMatching(gpr, higst, matches, stats);
      //higst = null;
      
      matches.clear();
      stats.println("\n***** keyword tree:" + msFile);
      ExactMatchingOld.exactMatching(gpr, sequences, matches);
      
    }
    stats.close();
  }
  

  /**
   * Mutation benchmark
   * @param p The parameter object
   * @param db The path to the database
   * @param statsFile the path to the statistic files.
   */
  public static void controlledMutationSearch(Parameters p, String db, String statsFile) {
    
    ProteinFastaSequence sequences = new ProteinFastaSequence(db);
    int notFound = 0;
    PrintWriter stats = null;
    try {
      stats = new PrintWriter(statsFile);
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-9);
    }
    
    // assume the input is an MS2 file for now
    for (String msFile : p.specFiles()) {
      
      String resultFile = p.getOutFileName() + ".grc";
      System.out.println("Processing " + resultFile);
      
      GappedPeptideResults gpr = new MSGDResultFileParser(resultFile).iterator().next();
      
      // maps the sequence id to the list of matches
      HashMap<Integer,ArrayList<MatchObject>> matches = new HashMap<Integer,ArrayList<MatchObject>>();
      
      stats.println("\n***** keyword tree:" + msFile);
      ExactMatchingOld.exactMatching(gpr, sequences, matches);
      
      SpectraIterator it = null;
      try {
        it = new SpectraIterator(msFile, new MS2SpectrumParser());
      }
      catch (IOException ioe) {
        System.err.println("Invalid spectrum file: " + msFile);
        System.err.println(ioe);
        System.exit(-1);
      }
      
      ArrayList<SpectrumMatches> specMatches = new ArrayList<SpectrumMatches>();
      int specCount = 0;
      // maps the query id to the list of exact match objects
      HashMap<Integer,ArrayList<MatchObject>> queryMatches = new HashMap<Integer,ArrayList<MatchObject>>();
      // populate the specMatches ArrayList
      while (it.hasNext()) {
        specCount++;
        Spectrum spec = it.next();
        queryMatches.clear();
        
        //System.out.println("Matches for spectrum " + spec.getScanNum());
        if (matches.containsKey(spec.getScanNum())) {
          // there are matches for this spectrum
          queryMatches.clear();
          /*
          SpectrumMatches sm = new SpectrumMatches(spec);
          //System.out.println(sm.getScoredSpec() + " Charge " + spec.getCharge());
          
          for (MatchObject emo : matches.get(spec.getScanNum())) {
            if (!queryMatches.containsKey(emo.getQueryIndex())) {
              queryMatches.put(emo.getQueryIndex(), new ArrayList<MatchObject>());
            }
            queryMatches.get(emo.getQueryIndex()).add(emo);
          }
          
          for (int queryIndex : queryMatches.keySet()) {
            sm.addMatches(queryMatches.get(queryIndex));
          }
          specMatches.add(sm); */
        }
        else {
          notFound++;
          //System.out.println("No matches"); 
        }
        //System.out.println();
      }
      System.out.println("Total spectra " + specCount + ". " + notFound + " with no matches.");
      

      // select spectrum match objects that do not overlap
      System.out.println("Total spectrum match objects before filtering " + specMatches.size());
      ArrayList<SpectrumMatches> filteredMatches = selectNonOverlappingMatches(specMatches);
      System.out.println("Total filtered (non-overlapping) spectrum matches objects " + filteredMatches.size());
      
      // Map the coordinates to exact match objects. This will eliminate redundant matches
      HashMap<Long,ExactMatchObject> coors2ExactMatches = new HashMap<Long,ExactMatchObject>();
      //HashMap<Long,Integer> coors2SpecId = new HashMap<Long,Integer>();

      for (SpectrumMatches sm : filteredMatches) {
        //System.out.println(sm.getScanNumber() + " spectrum break " + sm.getSequencesCount());
        for (MatchObject mo : sm.getMatches()) {
          System.err.println("Error tryign to get coordinates of match object " + mo.getStringWithCoordinates());
          //coors2ExactMatches.put(((ExactMatchObject)mo).getCoors(), (ExactMatchObject)mo);
          //coors2SpecId.put(((ExactMatchObject)mo).getCoors(), sm.getScanNumber());
        }
        //System.out.println(sm);
      }
      
      
      

      
      // mutate the database, randomly
      sequences.makeModifiable();
      HashMap<Long,String> coors2Mutations = new HashMap<Long,String>();
      HashMap<Long,Long> coors2Positions = new HashMap<Long,Long>();
      generateRandomPositions(sequences, filteredMatches, coors2Mutations, coors2Positions);
      

  
      // create the new gpr
      GappedPeptideResults fgpr = new GappedPeptideResults();
      /*
      for (SpectrumMatches sm : filteredMatches) {
        /*
        Spectrum s = sm.getSpectrum();
        fgpr.addSpectrum(s.getScanNum(), s.getScanNum(), msFile, s.getPrecursorPeak().getMz(), s.getCharge(), "null");
        for (int sequenceIndex=0; sequenceIndex < sm.getSequencesCount(); sequenceIndex++) {
          for (MatchObject mo : sm.getMatchesForSequenceAt(sequenceIndex)) {
            fgpr.addSequence(s.getScanNum(), mo.getQuery());
          }
        }
      }
      */
      
      // use the mutation tolerant search
      HashMap<Integer,ArrayList<MutMatchObject>> mutatedMatches = new HashMap<Integer,ArrayList<MutMatchObject>>();
      //oneMutationMatching(fgpr, sequences, mutatedMatches, "/home/jung/Desktop/mutationMatches.txt");
      
      int nonUniqueMatches = 0;
      TreeMap<Long,MutMatchObject> uniqueMatches = new TreeMap<Long,MutMatchObject>();
      for (int specId : mutatedMatches.keySet()) {
        for (MutMatchObject smmo : mutatedMatches.get(specId)) {
          long coors = 0;//smmo.getCoors();
          if (!uniqueMatches.containsKey(coors)) {
            uniqueMatches.put(coors, smmo);
          }
          nonUniqueMatches++;
        }
      }
      
      int nonRecoveredMatches = 0;
      
      /*
      int incorrectMatches = 0;
      for (SingleMutationMatchObject smmo : uniqueMatches.values()) {
        if (!Peptide.isCorrect(smmo.getMatchAsString(), smmo.getQuery())) {
          incorrectMatches++;
        }
        //System.out.printf("%s %s\n", smmo.getStringWithCoordinates(), smmo.toString());
      }
      */
      
      
      // verification
      System.out.println();
      for (long keyCoor : coors2ExactMatches.keySet()) {
        if (!uniqueMatches.containsKey(keyCoor)) {
          System.out.println(coors2ExactMatches.get(keyCoor).getStringWithCoordinates() + " not found");
          System.out.println(coors2Positions.get(keyCoor) + " ; " + coors2Mutations.get(keyCoor));
          System.out.println(coors2ExactMatches.get(keyCoor).getQuery() + " - " + coors2ExactMatches.get(keyCoor).getQueryIndex());
          //System.out.println("Found matches:");
          
          //for (SingleMutationMatchObject smmo : mutatedMatches.get(coors2SpecId.get(keyCoor))) {
          //  System.out.println(smmo.getQuery() + " - " +smmo.getStringWithCoordinates() + " - " + smmo.toString());
          //}
          
          //System.out.println();
          //System.exit(-9);
          nonRecoveredMatches++;
        }
      }
      
      System.out.println();
      System.out.printf("Total spectra queried %d\n", fgpr.getSpectrumCount());
      System.out.printf("Total sequences queried %d\n", fgpr.getSequenceCount());
      System.out.printf("Total non-unique matches found %d\n", nonUniqueMatches);
      System.out.printf("Total unique matches found %d\n", uniqueMatches.size());
      System.out.printf("Total not recovered matches %d\n", nonRecoveredMatches);
      //System.out.printf("Total incorrect matches %d\n", incorrectMatches);
    }
    stats.close();
  }
  
  
  private static void generateRandomPositions(ProteinFastaSequence sequences,
                                              ArrayList<SpectrumMatches> matches,
                                              HashMap<Long,String> coors2Mutations,
                                              HashMap<Long,Long> coors2Positions) {
    Random r = new Random();
    for (SpectrumMatches sm : matches) {
      for (MatchObject mo : sm.getMatches()) {
        long start = mo.getStart(), end = mo.getEnd();
        long position = start+r.nextInt((int)(end-start));
        //long position = start+3;
        
        char og = sequences.getCharAt(position);
        char mutation = 'x';
        //int aaIndex = 0;
        do {
          //mutation = AminoAcid.getStandardAminoAcids()[aaIndex++].getResidue();
          mutation = AminoAcid.getStandardAminoAcids()[r.nextInt(20)].getResidue();
        } while (AminoAcid.getStandardAminoAcid(og).getNominalMass()==AminoAcid.getStandardAminoAcid(mutation).getNominalMass());
        //System.out.printf("Mutating from %c to %c at %d\n", og, mutation, position);
        sequences.set(position, mutation);
        
        System.err.println("Error tryign to get coordinates of match object " + mo.getStringWithCoordinates());
        //coors2Mutations.put(((ExactMatchObject)mo).getCoors(), String.format("%c -> %c", og, mutation));
        //coors2Positions.put(((ExactMatchObject)mo).getCoors(), position-start);
        
      }
    }
  }

  /*
  private static boolean elegibleQuery(ArrayList<Integer> query) {
    int totalMass = 0;
    for (int mass : query) totalMass += mass;
    
    int cumMass = 0, threshold = Constants.MIN_QUERY_MASS;
    for (int mass : query) {
      // this will make the query fail if mutation falls here 
      if (cumMass < threshold && totalMass-cumMass-mass < threshold) return false;
    }
    
    return true;
  }
  */
  
  private static ArrayList<SpectrumMatches> selectNonOverlappingMatches(ArrayList<SpectrumMatches> matches) {
    
    // start and end positions 
    TreeMap<Long,Long> seenPositions = new TreeMap<Long,Long>();
    
    ArrayList<SpectrumMatches> results = new ArrayList<SpectrumMatches>();
    
    // do a simple first come, first reserve proteome position approach
    for (SpectrumMatches sm : matches) {
      boolean approved = true;
      for (MatchObject mo : sm.getMatches()) {
        long start = mo.getStart(), end = mo.getEnd();
        
        // check that the start index doesn't land inside a segment
        Long lower = seenPositions.floorKey(start);
        if (lower!=null) {
          if (seenPositions.get(lower) > start) {
            approved = false;
            break;
          }
        }
        
        // check that the end index doesn't land inside a segment
        Long upper = seenPositions.floorKey(end-1);
        if (upper!=null) {
          if (seenPositions.get(upper) >= end) {
            approved = false;
            break;
          }
        }

        /*
        if (!elegibleQuery(mo.getQuery())) {
          approved = false;
          break;
        }*/
        
        if (!approved) break;
      }
      
      if (approved) {
        results.add(sm);
        for (MatchObject mo : sm.getMatches()) {
          System.err.println("Broken method... cannot store " + mo.toString());
          //seenPositions.put(mo.getStart(), mo.getEnd());
        }
      }
    }
    
    return results;
  }
  
  
  
  public static void simpleMatch(Param p) {
    
    // try to process all the files ending with grc in the directory
    File dirObj = new File(p.input);
    for (String filename : dirObj.list()) {
      if (filename.endsWith(".grc")) {
        runSearch(p.input + "/" + filename, p.db, p.output, p.mutatedOutput);
        //break;
      }
    }
  }
  
  
  
  public static void main(String[] args) {
    
    String userHome = System.getProperty("user.home");
    
    String mutatedResultFile = userHome + "/Desktop/mutatedMatches.txt";
    String exactResultFile = userHome + "/Desktop/exactMatches.txt";
    
    //String dir = userHome + "/Desktop/liverResults";
    //String dir = userHome + "/Desktop/miniLiverResults";
    //String dir = userHome + "/Data/Results/HumanLiver";
    //String input = userHome + "/Data/Spectra/SOneSpectra/Yufeng_LTQ_15um_4000psi_50ngShew_dd_4_dta.ms2";
    String input = userHome + "/Data/Spectra/SOneSpectra/Shew054-160_22aug01_saturn_c2_400-2000_dta.ms2";
    
    //String humanProteinDB = userHome + "/Data/Databases/ensembl/ensemblCRCh37proteins.fasta";
    //String db = userHome + "/Data/Databases/human/ipi.HUMAN.v3.72.fasta";
    //String humanProteinDB = userHome+"/Data/Databases/yeast_nr050706.fasta";
    String db = userHome + "/Data/Databases/ShewDB/SOne_uniprot_plus_contaminants.fasta";

    //String msgdParams = userHome + "/Software/msgd/SOne_tiny.params";
    String msgdParams = userHome + "/Software/msgd/all.params";
    
    //String statsFile = userHome + "/Desktop/stats.txt";
    
    Parameters p = new Parameters(msgdParams, true);
    long time = System.currentTimeMillis();
    NewMSGappedDictionary.run(p, true);
    System.out.println("Time spent running MSGD " + (System.currentTimeMillis()-time) / 1000.0);
    
    Param param = new Param();
    param.input = input;
    param.output = exactResultFile;
    param.mutatedOutput = mutatedResultFile;
    param.db = db;
    
    //System.out.println(p.outFileName());
    //for (File specFile : p.specFiles())
    //System.out.println(specFile);

    //simpleMatch(p);
    
    //compareKeywordSuffixTree(p, db, statsFile);
    
    //controlledMutationSearch(p, db, statsFile);
  }
  
}
