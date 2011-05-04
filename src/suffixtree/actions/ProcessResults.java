package suffixtree.actions;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;

import cyclic.Cluster1D;
import cyclic.Point1D;

import msdictionary.Codon;
import msutil.AminoAcid;

import sequences.FastaSequence;
import sequences.ProteinFastaSequences;
import sequences.Sequence;
import suffixtree.Constants;
import suffixtree.results.ModdedResult;

/**
 * This class contains the methods to post-process the results
 * @author jung
 *
 */
public class ProcessResults {

  
  private static class Result implements Comparable<Result> {
    private long start;
    private long end;
    private long position;
    private char original;
    private char mutation;
    private float offset;
    private float prob;
    private boolean valid;
    private String annotation;
    private String id;
    private String peptide;
    private String codon;
    private String filename;
    private int scanNum;
    private int score;
    
    private Result(String line) {
      
      String[] tokens = line.split("\t");
      
      //System.out.println(line);
      
      // parse the fields
      String[] filenameTokens = tokens[0].split("/");
      int scanNum = Integer.parseInt(tokens[1]);
      String id = tokens[5];
      String annotation = tokens[6];
      float prob = Float.parseFloat(tokens[8]);
      String peptide = tokens[9];
      long start = Long.parseLong(tokens[10]);
      long end = Long.parseLong(tokens[11]);
      long position = Long.parseLong(tokens[12]);
      char original = tokens[13].charAt(0);
      char mutation = tokens[14].charAt(0);
      float offset = Float.parseFloat(tokens[15]);
      
      //this.line = line;
      this.score = Integer.parseInt(tokens[7]);
      this.scanNum = scanNum;
      this.start = start;
      this.end = end;
      this.original = original;
      this.mutation = mutation;
      this.offset = offset;
      this.position = position;
      this.prob = prob;
      this.valid = true;
      this.id = id;
      this.annotation = annotation;
      this.codon = "---";
      this.peptide = peptide;
      this.filename = filenameTokens[filenameTokens.length-1];
    }
    
    @Override
    public String toString() {
      return String.format("%s\t%d\t%d\t%d\t%s\t%d\t%.2e\t%c\t%c\t%f\t%f\t%s\t%s\t%s", filename, scanNum, start, end, annotation, score, prob, original, mutation, offset, delta(), codon, id, peptide);
    }

    @Override
    public int compareTo(Result o) {
      // reverse order sort. smaller first
      if (this.start > o.start) return 1;
      if (o.start > this.start) return -1;
      return 0;
    }
    
    private boolean isInsertion() {
      return this.original==suffixtree.Constants.EMPTY_AA;
    }
    
    private int getOffset() {
      return Integer.parseInt(this.annotation.split("\\s")[3]);
    }
    
    private int getShift() {
      return Integer.parseInt(this.annotation.split("\\s")[1]);
    }
    
    private float getMassError() { return this.offset; }
    
    private String getName() {
      return this.annotation.split("\\s")[0];  
    }
    
    private boolean isReversed() {
      return Integer.parseInt(this.annotation.split("\\s")[2])==1;
    }
    
    private float delta() {
      float oMass = 0.0f, mMass = 0.0f;
      if (original!=suffixtree.Constants.EMPTY_AA) oMass = Constants.AA.getAminoAcid(original).getMass();
      if (mutation!=suffixtree.Constants.EMPTY_AA) mMass = Constants.AA.getAminoAcid(mutation).getMass();
      return offset + mMass - oMass;
      //return offset + oMass - mMass;
    }
  }
  
  
  private static void printStatistics(PrintWriter out, ArrayList<Result> results) {
    
    /**
     * Class to make things sortable in reverse order
     * @author jung
     *
     */
    class Pair implements Comparable<Pair>{
      private String s;
      private float f;
      
      private Pair(String s, float f) {
        this.s = s;
        this.f = f;
      }

      @Override
      public int compareTo(Pair o) {
        if (this.f > o.f) return -1;
        if (o.f > this.f) return 1;
        return 0;
      }
      
      @Override
      public String toString() {
        return String.format("%s\t%.2f%%", this.s, this.f);
      }
    }
    
    final int MAX_SHIFT = 3;
    
    // collect the codon statistics
    out.print("\nCodons distribution\n");
    HashMap<String,Integer> codons = new HashMap<String,Integer>();
    for (Result r : results) {
      if (!codons.containsKey(r.codon)) {
        codons.put(r.codon, 1);
      }
      else {
        codons.put(r.codon, codons.get(r.codon)+1);
      }
    }
    ArrayList<Pair> items = new ArrayList<Pair>(); 
    for (String codon : codons.keySet()) {
      items.add(new Pair(codon, codons.get(codon) * 100.0f / results.size()));
    }
    Collections.sort(items);
    for (Pair item : items) {
      out.println(item);
    }
    
    // collect the neighboring amino acid statistics
    out.printf("\nNeighboring amino acids up to +/-%d positions\n", MAX_SHIFT);
    HashMap<Character,Integer> vecinos = new HashMap<Character,Integer>();
    int totalCount = 0;
    for (Result r : results) {
      //System.out.print(r.id + "\t");
      int offset = (int)(r.position - r.start);
      for (int shift=-MAX_SHIFT; shift<0; shift++) {
        int index = offset+shift;
        if (index<0) continue;
        char aa = r.id.charAt(index);
        if (!vecinos.containsKey(aa)) {
          vecinos.put(aa, 1);
        }
        else {
          vecinos.put(aa, vecinos.get(aa)+1);
        }
        totalCount++;
      }
      for (int shift=0; shift<MAX_SHIFT; shift++) {
        int index = offset+shift+6;
        if (index>=r.id.length()) break;
        char aa = r.id.charAt(index);
        //System.out.print(aa + "\t");
        if (!vecinos.containsKey(aa)) {
          vecinos.put(aa, 1);
        }
        else {
          vecinos.put(aa, vecinos.get(aa)+1);
        }
        totalCount++;
      }
      //System.out.println();
    }
    items.clear(); 
    for (char aa : vecinos.keySet()) {
      items.add(new Pair(aa+"", vecinos.get(aa) * 100.0f / totalCount));
    }
    Collections.sort(items);
    for (Pair item : items) {
      out.println(item);
    }
    
  }
  
  
  
  private static void makeMutationTable(String dir, ArrayList<Result> results) {
    HashMap<Integer,Character> index2aa = new HashMap<Integer,Character>();
    int aaIndex = 1;
    index2aa.put(0, suffixtree.Constants.EMPTY_AA);
    Iterator<AminoAcid> it = Constants.AA.iterator();
    while (it.hasNext()) {
      AminoAcid aa = it.next();      
      index2aa.put(aaIndex++, aa.getResidue());
    }
    
    HashMap<String,ArrayList<Result>> entries = new HashMap<String,ArrayList<Result>>();
    for (Result result : results) {
      String key = result.original + "" + result.mutation;
      if (!entries.containsKey(key)) {
        entries.put(key, new ArrayList<Result>());
      }
      entries.get(key).add(result);
    }

    String filename = dir + "/table.html";
    if (!new File(dir).exists()) {
      // create the directory
      boolean result = new File(dir).mkdir();
      if (!result) System.out.println("Failed to create directory");
    }
    
    // the running totals for each column
    int[] totals = new int[index2aa.size()];
    
    StringBuffer sb = new StringBuffer();
    sb.append("<table border='1'><tr><th></th>");
    // header
    sb.append("<th>*(0)</th>");
    for (int j=1; j<index2aa.size(); j++) {
      AminoAcid aa = Constants.AA.getAminoAcid(index2aa.get(j));
      sb.append(String.format("<th>%c(%d)</th>", index2aa.get(j), aa.getNominalMass()));
    }
    sb.append("<th>Total</th></tr>\n");
    
    // table body
    for (int i=0; i<index2aa.size(); i++) {
      sb.append("<tr>");
      if (i==0) {
        sb.append("<th>*(0)</th>");
      }
      else {
        AminoAcid aa = Constants.AA.getAminoAcid(index2aa.get(i));
        sb.append(String.format("<th>%c(%d)</th>", index2aa.get(i), aa.getNominalMass()));
      }
      int total = 0;
      for (int j=0; j<index2aa.size(); j++) {
        String key = index2aa.get(i) + "" + index2aa.get(j);
        if (entries.containsKey(key)) {
          
          // create the entry
          String subFilename = dir+"/cell"+key.replace('*', '_')+".txt";
          try {
            PrintWriter p = new PrintWriter(subFilename);
            Collections.sort(entries.get(key));
            for (Result r : entries.get(key)) {
              p.println(r);
            }
            printStatistics(p, entries.get(key));
            p.close();
          }
          catch (IOException ioe) {
            System.err.println(ioe);
          }
          total += entries.get(key).size();
          totals[j] += entries.get(key).size();
          sb.append(String.format("<td><a href='%s' target='_blank'>%d</a></td>", "cell"+key.replace('*', '_')+".txt", entries.get(key).size()));
        }
        else {
          sb.append("<td>0</td>");
        }
      }
      sb.append(String.format("<td>%d</td>", total));
      sb.append("</tr>\n");
    }
    
    // the totals column
    sb.append("<tr><th>Total</th>");
    int grandTotal = 0;
    for (int total : totals) {
      sb.append(String.format("<td>%d</td>", total));
      grandTotal += total;
    }
    sb.append(String.format("<td>%d</td></tr>", grandTotal));
    sb.append("</table>");
    
    try {
      PrintWriter out = new PrintWriter(filename);
      out.println("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">");
      out.println("<head><title>Mutation Table</title></head><body>");
      out.println(sb);
      out.println("</body></html>");
      out.close();
    }
    catch (IOException ioe) {
      System.err.println(ioe);
    }
    
  }
  
  /*
  private int translateCoordinate(int pos, int shift, int offset, boolean isReversed) {
    if (!isReversed) {
      return pos*3+offset-3;
    }
    return 0;  
  }
  */

  /**
   * Look the codon back from the DNA sequence
   * @param results the list of items to populate
   * @param dna the DNA sequence
   */
  private static void populateCodons(ArrayList<Result> results, Sequence dna) {
    
    for (Result result : results) {
      
      String entry = dna.getMatchingEntry(result.getName());
      int shift = result.getOffset();
      int offset = result.getShift();
      if (!result.isReversed()) {
        if (!result.isInsertion()) {
          result.codon = entry.substring((int)(result.position*3+offset-3), (int)(offset+result.position*3)).toUpperCase();
          if (Codon.translate(result.codon)!=result.original) {
            System.err.printf("Incorrect forward codons %c\t%s\t%s\t%d\n", result.original, result.codon, result.getName(), offset);
          }
        }
      }
      else {
        if(!result.isInsertion()) {
          // find the first non-ATCG position after the offset
          int end=offset;
          for (; end<entry.length(); end++) {
            char letter = Character.toUpperCase(entry.charAt(end));
            if (!(letter=='A' || letter=='C' || letter=='G' || letter=='T')) {
              break;
            }
          }
          long position = end - shift - 3*result.position;
          StringBuffer sb = new StringBuffer(entry.substring((int)position, (int)(position+3)).toUpperCase()).reverse();
          StringBuffer comp = new StringBuffer();
          for (int i=0; i<sb.length(); i++) {
            comp.append(Codon.complement(sb.charAt(i)));
          }
          result.codon = comp.toString();
          if (Codon.translate(result.codon)!=result.original) {
            System.err.printf("Incorrect reverse codons %c\t%s\t%d\t%d\t%s\n", result.original, result.codon, offset, shift, result.getName());
          }
        }
      }
    }
  }
  
  
  
  private static ArrayList<Result> keepByDelta(ArrayList<Result> results, float cutOff) {
    ArrayList<Result> retItems = new ArrayList<Result>();
    for (Result result : results) {
      if (Math.abs(result.getMassError()) <= cutOff) retItems.add(result);
    }
    return retItems;
  }
  
  
  
  /**
   * Filter out all those items that do not have the mutation at the beginning
   * @param results the list of results
   * @return the filtered list of results
   */
  private static ArrayList<Result> keepFirstPosition(ArrayList<Result> results) {
    ArrayList<Result> retItems = new ArrayList<Result>();
    for (Result result : results) {
      if (result.position==result.start) retItems.add(result);
    }
    return retItems;
  }
  
  
  private static ArrayList<Result> probFilter(ArrayList<Result> results, float prob) {
    ArrayList<Result> retItems = new ArrayList<Result>();
    for (Result result : results) {
      if (result.prob <= prob) retItems.add(result);
    }
    return retItems;
  }
  
  private static ArrayList<Result> scoreFilter(ArrayList<Result> results, int score) {
    ArrayList<Result> retItems = new ArrayList<Result>();
    for (Result result : results) {
      if (result.score >= score) retItems.add(result);
    }
    return retItems;
  }
  
  
  /**
   * Eliminate all those entries that are repeated. Mutations have preferences over
   * insertions, and lower probability is better. 
   * @param results the list of items to filter
   * @return the filtered items
   */
  private static ArrayList<Result> eliminateRedundant(ArrayList<Result> results) {
    
    // group the items by the start position
    HashMap<Long,ArrayList<Result>> groupedByStart = new HashMap<Long,ArrayList<Result>>();
    for (Result result : results) {
      if (!groupedByStart.containsKey(result.start)) {
        groupedByStart.put(result.start, new ArrayList<Result>());
      }
      groupedByStart.get(result.start).add(result);
    }
    
    for (long startPos : groupedByStart.keySet()) {
      HashMap<Character,Result> groupedByMutation = new HashMap<Character,Result>();
      // group by mutation
      for (Result result : groupedByStart.get(startPos)) {
        if (!groupedByMutation.containsKey(result.mutation)) {
          groupedByMutation.put(result.mutation, result);
        }
        else {
          Result previous = groupedByMutation.get(result.mutation);
          if (previous.isInsertion()) {
            if (result.isInsertion()) {
              if (result.prob < previous.prob) {
                groupedByMutation.put(result.mutation, result);
              }
            }
            else {
              groupedByMutation.put(result.mutation, result);
            }
          }
          else {
            if (!result.isInsertion()) {
              if (result.prob < previous.prob) {
                groupedByMutation.put(result.mutation, result);
              }
            }
          }
        }
      }
      
      // record, reset and restore selected items
      ArrayList<Boolean> selectedValues = new ArrayList<Boolean>();
      ArrayList<Result> selected = new ArrayList<Result>(groupedByMutation.values());
      for (Result result : selected) {
        selectedValues.add(result.valid);
      }
      for (Result result : groupedByStart.get(startPos)) {
        //System.out.println("Set false");
        result.valid = false;
      }
      for (int index=0; index<selectedValues.size(); index++) {
        selected.get(index).valid = selectedValues.get(index);
      }
    }
    
    
    // group the items by the end position
    HashMap<Long,ArrayList<Result>> groupedByEnd = new HashMap<Long,ArrayList<Result>>();
    for (Result result : results) {
      if (!groupedByEnd.containsKey(result.end)) {
        groupedByEnd.put(result.end, new ArrayList<Result>());
      }
      groupedByEnd.get(result.end).add(result);
    }
    
    for (long endPos : groupedByEnd.keySet()) {
      HashMap<Character,Result> groupedByMutation = new HashMap<Character,Result>();
      // group by mutation
      for (Result result : groupedByEnd.get(endPos)) {
        if (!groupedByMutation.containsKey(result.mutation)) {
          groupedByMutation.put(result.mutation, result);
        }
        else {
          Result previous = groupedByMutation.get(result.mutation);
          if (previous.isInsertion()) {
            if (result.isInsertion()) {
              if (result.prob < previous.prob) {
                groupedByMutation.put(result.mutation, result);
              }
            }
            else {
              groupedByMutation.put(result.mutation, result);
            }
          }
          else {
            if (!result.isInsertion()) {
              if (result.prob < previous.prob) {
                groupedByMutation.put(result.mutation, result);
              }
            }
          }
        }
      }
      
      // record, reset and restore selected items
      ArrayList<Boolean> selectedValues = new ArrayList<Boolean>();
      ArrayList<Result> selected = new ArrayList<Result>(groupedByMutation.values());
      for (Result result : selected) {
        selectedValues.add(result.valid);
      }
      for (Result result : groupedByEnd.get(endPos)) {
        result.valid = false;
      }
      for (int index=0; index<selectedValues.size(); index++) {
        selected.get(index).valid = selectedValues.get(index);
      }
    }
    
    ArrayList<Result> filteredResults = new ArrayList<Result>();
    for (Result current : results) {
      if (current.valid) filteredResults.add(current);
    }
    return filteredResults;
  }
  
  
  public static void processMutatedResults() {
    
    String userHome = System.getProperty("user.home");

    String resultFile = userHome + "/Data/Spectra/Ecoli/ORG045_LTQ_Orb_0/allResultsMutated.txt";
    //String translatedDbFile = userHome + "/Data/Databases/Ecoli/pro/translated.fasta";
    String dnaDbFile = userHome + "/Data/Databases/Ecoli/gen/E_coli_BL21.EB1.dna.toplevel.fasta";
    String outDir = userHome + "/Data/Spectra/Ecoli/ORG045_LTQ_Orb_0/allResultsMutationTable";

    Sequence dnaDb = new FastaSequence(dnaDbFile);
    
    ArrayList<Result> results = new ArrayList<Result>();
    
    try {
      BufferedReader in = new BufferedReader(new FileReader(resultFile));
      
      String line;
      while ((line = in.readLine())!=null) {
        Result current = new Result(line);
        results.add(current);
      }
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    int initialSize = results.size();
    ArrayList<Result> filteredResults = results;
    filteredResults = keepFirstPosition(filteredResults);
    filteredResults = keepByDelta(filteredResults, 0.01f); 
    filteredResults = eliminateRedundant(filteredResults);
    populateCodons(filteredResults, dnaDb);
    makeMutationTable(outDir, filteredResults);
    System.out.printf("Previous count: %d. Filtered count: %d\n", initialSize, filteredResults.size()); 
  }
  
  public static void processMultipleMutatedResultsAvar() {
    
    String userHome = System.getProperty("user.home");

    String[] resultFiles =  {userHome + "/Data/Spectra/Avar/ORG013_LTQ_0/allResultsMutated.txt"};
                           //  userHome + "/Data/Spectra/Avar/ORG013_LTQ_1/allResultsMutated.txt"};
    String dnaDbFile = userHome + "/Data/Databases/Avar/gen/Avar.fasta";
    String outDir = userHome + "/Data/Spectra/Avar/allResultsMutationTable";

    Sequence dnaDb = new FastaSequence(dnaDbFile);
    
    ArrayList<Result> results = new ArrayList<Result>();
    
    try {
      for (String resultFile : resultFiles) {
        BufferedReader in = new BufferedReader(new FileReader(resultFile));
        
        String line;
        while ((line = in.readLine())!=null) {
          Result current = new Result(line);
          results.add(current);
        }
      }
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    int initialSize = results.size();
    ArrayList<Result> filteredResults = results;
    filteredResults = probFilter(filteredResults, 1e-11f);
    filteredResults = keepFirstPosition(filteredResults);
    filteredResults = keepByDelta(filteredResults, 0.05f); 
    filteredResults = eliminateRedundant(filteredResults);
    populateCodons(filteredResults, dnaDb);
    makeMutationTable(outDir, filteredResults);
    System.out.printf("Previous count: %d. Filtered count: %d\n", initialSize, filteredResults.size());    
  }
  
  
  public static void processMultipleMutatedResultsCsp() {
    
    String userHome = System.getProperty("user.home");

    String[] resultFiles =  {userHome + "/Data/Spectra/Csp/ORG033_LTQ_Orb_0/resultsMutated6.txt"};
                             //userHome + "/Data/Spectra/Csp/ORG033_LTQ_0/allResultsMutated.txt"};
    String dnaDbFile = userHome + "/Data/Databases/Csp/gen/Csp.fasta";
    String outDir = userHome + "/Data/Spectra/Csp/mutationTable";

    Sequence dnaDb = new FastaSequence(dnaDbFile);
    
    ArrayList<Result> results = new ArrayList<Result>();
    
    try {
      for (String resultFile : resultFiles) {
        BufferedReader in = new BufferedReader(new FileReader(resultFile));
        
        String line;
        while ((line = in.readLine())!=null) {
          if (line.split("\\s").length<14) continue;
          Result current = new Result(line);
          results.add(current);
        }
      }
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    int initialSize = results.size();
    ArrayList<Result> filteredResults = results;
    filteredResults = probFilter(filteredResults, 1e-13f);
    filteredResults = keepFirstPosition(filteredResults);
    filteredResults = keepByDelta(filteredResults, 0.02f); 
    filteredResults = eliminateRedundant(filteredResults);
    populateCodons(filteredResults, dnaDb);
    makeMutationTable(outDir, filteredResults);
    System.out.printf("Previous count: %d. Filtered count: %d\n", initialSize, filteredResults.size());    
  }
  
  

  public static void processMultipleMutatedResultsAsp() {
    
    String userHome = System.getProperty("user.home");

    String[] resultFiles = {userHome + "/Data/Spectra/Asp/ORG014_LTQ_FT_0/resultsMutated.txt"};
                            //userHome + "/Data/Spectra/Asp/ORG014_LTQ_3/allResultsMutated.txt",
                            //userHome + "/Data/Spectra/Asp/ORG014_LTQ_2/allResultsMutated.txt",
                            //userHome + "/Data/Spectra/Asp/ORG014_LTQ_1/allResultsMutated.txt",
                            //userHome + "/Data/Spectra/Asp/ORG014_LTQ_0/allResultsMutated.txt"};
    String dnaDbFile = userHome + "/Data/Databases/Asp/gen/Asp.fasta";
    String outDir = userHome + "/Data/Spectra/Asp/mutationTable";

    Sequence dnaDb = new FastaSequence(dnaDbFile);
    
    ArrayList<Result> results = new ArrayList<Result>();
    
    // collect only the better score
    HashMap<String,Integer> scores = new HashMap<String,Integer>();
    
    try {
      for (String resultFile : resultFiles) {
        BufferedReader in = new BufferedReader(new FileReader(resultFile));
        
        String line;
        while ((line = in.readLine())!=null) {
          Result current = new Result(line);
          String key = current.filename+"$$"+current.scanNum;
          if (scores.containsKey(key)) {
            if (current.score > scores.get(key)) {
              scores.put(key, current.score);
            }
          }
          else {
            scores.put(key, current.score);
          }
          results.add(current);
        }
      }
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    // group the scores and print them in order
    TreeMap<Integer,Integer> scoreHist = new TreeMap<Integer,Integer>();
    for (int score : scores.values()) {
      if (scoreHist.containsKey(score)) {
        scoreHist.put(score, scoreHist.get(score)+1);
      }
      else {
        scoreHist.put(score, 1);
      }
    }
    // print out the histogram
    for (int score : scoreHist.keySet()) {
      System.out.printf("%d\t%d\n", score, scoreHist.get(score));
    }
    
    
    int initialSize = results.size();
    ArrayList<Result> filteredResults = results;
    //filteredResults = probFilter(filteredResults, 1e-13f);
    filteredResults = scoreFilter(filteredResults, 159);
    filteredResults = keepFirstPosition(filteredResults);
    filteredResults = keepByDelta(filteredResults, 0.02f); 
    filteredResults = eliminateRedundant(filteredResults);
    populateCodons(filteredResults, dnaDb);
    makeMutationTable(outDir, filteredResults);
    System.out.printf("Previous count: %d. Filtered count: %d\n", initialSize, filteredResults.size());    
  }
  
  
  public static void collectOffsetDist(String[] dirs, int cutoff, String outdir) {
    HashMap<String,ModdedResult> results = new HashMap<String,ModdedResult>();
    try {
      for (String resultDir : dirs) {
        BufferedReader in = new BufferedReader(new FileReader(resultDir));
        String line;
        while ((line=in.readLine())!=null) {
          ModdedResult r = new ModdedResult(line);
          
          String peptide = r.getPeptide(); 
          if (peptide.charAt(0)!='K' && peptide.charAt(0)!='R') {
            //System.out.println(peptide);
            continue;
          }
          
          //if (r.getScore() >= cutoff) {
            String key = r.getFilepath() + r.getScanNumber();
            if (!results.containsKey(key)) {
              results.put(key, r);
            }
            else {
              //if (r.getScore() > results.get(key).getScore()) {
                // take the better score
                results.put(key, r);  
              //}
            }
          //}
        }
      }
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    // bin the offsets
    System.out.println("Number of items with passing score " + results.size());
    
    ArrayList<Point1D> points = new ArrayList<Point1D>();
    
    // make a histogram
    TreeMap<Integer,Integer> offsetCounts = new TreeMap<Integer,Integer>();
    int resultIndex = 0;
    ArrayList<ModdedResult> resultArray = new ArrayList<ModdedResult>();
    for (ModdedResult r : results.values()) {
      points.add(new Point1D(r.getDelta(), 1, resultIndex++));
      resultArray.add(r);
      int intOffset = r.getIntegerOffset();
      if (offsetCounts.containsKey(intOffset)) {
        offsetCounts.put(intOffset, offsetCounts.get(intOffset)+1);
      }
      else {
        offsetCounts.put(intOffset, 1);
      }
    }
    
    Set<Cluster1D> clusters = Cluster1D.cluster(points, 0.05f, 0);
    
    // print out the results / per offset
    try {
      PrintWriter out = new PrintWriter(outdir+"/offsetDist.html");
      out.println("<table border='1'><tr><th>Offset</th><th>Spectra</th><th>Peptides</th></tr>");
      for (Cluster1D c : clusters) {
        HashSet<String> sites = new HashSet<String>();
        // print out to a separate file
        String filename = String.format("Mods/file_%.3f.html", c.getCenter());
        String filepath = String.format("/home/jung/Data/Spectra/Sone/%s", filename);
        PrintWriter pw = new PrintWriter(filepath);
        pw.println("<table border='1'><tr><th>Peptide</th><th>Delta</th><th>Charge</th><th>Filename</th><th>Scan</th><th>Score</th><th>Prob</th></tr>");
        for (Point1D p : c.getPoints()) {
          ModdedResult r = resultArray.get(p.getIndex());
          sites.add(r.getModificationPositionKey());
          //pw.printf("<tr><td>%s</td><td>%.3f</td><td>%d</td><td>%s</td><td>%d</td><td>%d</td><td>%.3e</td></tr>", r.getPeptide(), r.getDelta(), r.getCharge(), r.getFilename(), r.getScanNumber(), r.getScore(), r.getProb());
          //pw.println(r.toString());
        }
        pw.println("</table>");
        pw.close();
        if (sites.size() >= 5)
        out.printf("<tr><td><a href='%s'>%.3f</a></td><td>%d</td><td>%d</td></tr>", filename, c.getCenter(), (int)c.getWeight(), sites.size());
      }
      out.println("</table>");
      
      /*
      for (int offset : offsetCounts.keySet()) {
        out.printf("%d\t%d\n", offset, offsetCounts.get(offset));
      }*/
      out.close();
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
  }

  
  
  public static void collectScoreDist(String[] dirs, String outfile) {
    HashMap<String,Integer> scores = new HashMap<String,Integer>();
    try {
      for (String resultDir : dirs) {
        BufferedReader in = new BufferedReader(new FileReader(resultDir));
        String line;
        while ((line=in.readLine())!=null) {
          ModdedResult r = new ModdedResult(line);
          String key = r.getFilepath() + r.getScanNumber();
          if (!scores.containsKey(key)) {
            //scores.put(key, r.getScore());
          }
          else {
            //if (r.getScore() > scores.get(key)) {
              // take the better score
            //  scores.put(key, r.getScore());  
            //}
          }
        }
      }
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }

    System.out.println("Number of items with scores " + scores.size());
    
    // make a histogram
    TreeMap<Integer,Integer> sortedScores = new TreeMap<Integer,Integer>();
    for (int score : scores.values()) {
      if (sortedScores.containsKey(score)) {
        sortedScores.put(score, sortedScores.get(score)+1);
      }
      else {
        sortedScores.put(score, 1);
      }
    }
    
    try {
      PrintWriter out = new PrintWriter(outfile);
      for (int score : sortedScores.keySet()) {
        out.printf("%d\t%d\n", score, sortedScores.get(score));
      }
      out.close();
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
  }

  
  public static void processMultipleMutatedResultsEcoli() {
    
    String userHome = System.getProperty("user.home");

    String[] resultFiles = {userHome + "/Data/Spectra/Ecoli/ORG045_LTQ_Orb_0/resultsMutated6.txt",
                            userHome + "/Data/Spectra/Ecoli/ORG045_LTQ_Orb_1/resultsMutated6.txt",
                            userHome + "/Data/Spectra/Ecoli/ORG045_LTQ_Orb_2/resultsMutated6.txt",
                            userHome + "/Data/Spectra/Ecoli/ORG045_LTQ_Orb_3/resultsMutated6.txt"};
    String dnaDbFile = userHome + "/Data/Databases/Ecoli/gen/E_coli_BL21.EB1.dna.toplevel.fasta";
    String outDir = userHome + "/Data/Spectra/Ecoli/MutationTable";

    Sequence dnaDb = new FastaSequence(dnaDbFile);
    
    ArrayList<Result> results = new ArrayList<Result>();
    
    try {
      for (String resultFile : resultFiles) {
        BufferedReader in = new BufferedReader(new FileReader(resultFile));
        
        String line;
        while ((line = in.readLine())!=null) {
          Result current = new Result(line);
          results.add(current);
        }
      }
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    int initialSize = results.size();
    ArrayList<Result> filteredResults = results;
    filteredResults = probFilter(filteredResults, 1e-13f);
    filteredResults = keepFirstPosition(filteredResults);
    filteredResults = keepByDelta(filteredResults, 0.02f); 
    filteredResults = eliminateRedundant(filteredResults);
    populateCodons(filteredResults, dnaDb);
    makeMutationTable(outDir, filteredResults);
    System.out.printf("Previous count: %d. Filtered count: %d\n", initialSize, filteredResults.size());    
  }
  
  
  
  public static void processMultipleModdedResultsSone() {
    String userHome = System.getProperty("user.home");
    
    String[] resultFiles = {userHome+"/Data/Spectra/Sone/LTQFT0/output6Modded.txt",
                            userHome+"/Data/Spectra/Sone/LTQFT1/output6Modded.txt",
                            userHome+"/Data/Spectra/Sone/LTQFT2/output6Modded.txt",
                            userHome+"/Data/Spectra/Sone/LTQFT3/output6Modded.txt",
                            userHome+"/Data/Spectra/Sone/LTQFT4/output6Modded.txt",
                            userHome+"/Data/Spectra/Sone/LTQFT5/output6Modded.txt"};
    String output = userHome+"/Data/Spectra/Sone/scoreDist.txt";
    collectScoreDist(resultFiles, output);
    
    // the reverse score distributions
    String[] resultFilesR = {userHome+"/Data/Spectra/Sone/LTQFT0/output6ModdedR.txt",
                             userHome+"/Data/Spectra/Sone/LTQFT1/output6ModdedR.txt",
                             userHome+"/Data/Spectra/Sone/LTQFT2/output6ModdedR.txt",
                             userHome+"/Data/Spectra/Sone/LTQFT3/output6ModdedR.txt",
                             userHome+"/Data/Spectra/Sone/LTQFT4/output6ModdedR.txt",
                             userHome+"/Data/Spectra/Sone/LTQFT5/output6ModdedR.txt"};
    String outputR = userHome+"/Data/Spectra/Sone/scoreDistR.txt";
    collectScoreDist(resultFilesR, outputR);
    
    // print out the offset counts given the cutoff
    String offsetDir = userHome+"/Data/Spectra/Sone";
    collectOffsetDist(resultFiles, 170, offsetDir);
  }
  
  
  
  public static void processMultipleMutatedResultsScervMito() {
    
    String userHome = System.getProperty("user.home");

    String[] resultFiles = {userHome + "/Data/Spectra/Scerv/ORG105_LTQ_Orb_0/resultsMitoMutated.txt"};
    String dnaDbFile = userHome + "/Data/Databases/Scerv/gen/chrmt.fasta";
    String outDir = userHome + "/Data/Spectra/Scerv/mutationTableMito";

    Sequence dnaDb = null;
    if (!new File(dnaDbFile).isDirectory()) {
      dnaDb = new FastaSequence(dnaDbFile);
    }
    else {
      dnaDb = new ProteinFastaSequences(dnaDbFile, true);
    }
    
    ArrayList<Result> results = new ArrayList<Result>();
    
    try {
      for (String resultFile : resultFiles) {
        BufferedReader in = new BufferedReader(new FileReader(resultFile));
        
        String line;
        while ((line = in.readLine())!=null) {
          if (line.split("\t").length < 12 ) continue;
          Result current = new Result(line);
          results.add(current);
        }
      }
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    int initialSize = results.size();
    ArrayList<Result> filteredResults = results;
    filteredResults = probFilter(filteredResults, 1e-11f);
    //filteredResults = keepFirstPosition(filteredResults);
    filteredResults = keepByDelta(filteredResults, 0.02f); 
    filteredResults = eliminateRedundant(filteredResults);
    populateCodons(filteredResults, dnaDb);
    makeMutationTable(outDir, filteredResults);
    System.out.printf("Previous count: %d. Filtered count: %d\n", initialSize, filteredResults.size());    
  }

  public static void processMultipleMutatedResultsScerv() {
    
    String userHome = System.getProperty("user.home");

    String[] resultFiles = {userHome + "/Data/Spectra/Scerv/ORG105_LTQ_Orb_0/resultsMutated6.txt",
                            userHome + "/Data/Spectra/Scerv/ORG105_LTQ_FT_0/resultsMutated6.txt"};
    String dnaDbFile = userHome + "/Data/Databases/Scerv/gen";
    String outDir = userHome + "/Data/Spectra/Scerv/mutationTable";

    Sequence dnaDb = null;
    if (!new File(dnaDbFile).isDirectory()) {
      dnaDb = new FastaSequence(dnaDbFile);
    }
    else {
      dnaDb = new ProteinFastaSequences(dnaDbFile, true);
    }
    
    ArrayList<Result> results = new ArrayList<Result>();
    
    try {
      for (String resultFile : resultFiles) {
        BufferedReader in = new BufferedReader(new FileReader(resultFile));
        
        String line;
        while ((line = in.readLine())!=null) {
          if (line.split("\t").length < 12 ) continue;
          Result current = new Result(line);
          results.add(current);
        }
      }
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    int initialSize = results.size();
    ArrayList<Result> filteredResults = results;
    filteredResults = probFilter(filteredResults, 1e-13f);
    //filteredResults = keepFirstPosition(filteredResults);
    filteredResults = keepByDelta(filteredResults, 0.02f); 
    filteredResults = eliminateRedundant(filteredResults);
    populateCodons(filteredResults, dnaDb);
    makeMutationTable(outDir, filteredResults);
    System.out.printf("Previous count: %d. Filtered count: %d\n", initialSize, filteredResults.size());    
  }


  public static void main(String[] args) {
    //processMutatedResults();
    
    //System.out.println("Mass of C is " + AminoAcid.getStandardAminoAcid('C').getMass());
    
    processMultipleMutatedResultsAsp();
    //processMultipleMutatedResultsEcoli();
    //processMultipleMutatedResultsCsp();
    //processMultipleMutatedResultsAvar();
    //processMultipleMutatedResultsScerv();
    //processMultipleMutatedResultsScervMito();
    //processMultipleModdedResultsSone();
  }
}
