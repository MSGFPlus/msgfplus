package suffixtree.actions;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;


import msdictionary.Codon;
import msutil.AminoAcid;

import sequences.FastaSequence;
import sequences.ProteinFastaSequences;
import sequences.Sequence;
import suffixtree.Constants;
import suffixtree.results.MutatedResult;





/**
 * This class collects the mutations and group them accordingly given a result
 * file
 * @author jung
 *
 */
public class CollectMutations {

  
  private static void printStatistics(PrintWriter out, ArrayList<MutatedResult> results) {
    
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
    for (MutatedResult r : results) {
      if (!codons.containsKey(r.getCodon())) {
        codons.put(r.getCodon(), 1);
      }
      else {
        codons.put(r.getCodon(), codons.get(r.getCodon())+1);
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
    for (MutatedResult r : results) {
      //System.out.print(r.id + "\t");
      int offset = (int)(r.getMutationPosition() - r.getStartPosition());
      for (int shift=-MAX_SHIFT; shift<0; shift++) {
        int index = offset+shift;
        if (index<0) continue;
        char aa = r.getPeptide().charAt(index);
        if (!vecinos.containsKey(aa)) {
          vecinos.put(aa, 1);
        }
        else {
          vecinos.put(aa, vecinos.get(aa)+1);
        }
        totalCount++;
      }
      for (int shift=0; shift<MAX_SHIFT; shift++) {
        int index = offset+shift+1;
        if (index>=r.getPeptide().length()) break;
        char aa = r.getPeptide().charAt(index);
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
  
  
  
  private static void makeMutationTable(String dir, ArrayList<MutatedResult> results) {
    HashMap<Integer,Character> index2aa = new HashMap<Integer,Character>();
    int aaIndex = 1;
    index2aa.put(0, suffixtree.Constants.EMPTY_AA);
    Iterator<AminoAcid> it = Constants.AA.iterator();
    while (it.hasNext()) {
      AminoAcid aa = it.next();      
      index2aa.put(aaIndex++, aa.getResidue());
    }
    
    HashMap<String,ArrayList<MutatedResult>> entries = new HashMap<String,ArrayList<MutatedResult>>();
    for (MutatedResult result : results) {
      String key = result.getOriginalAA() + "" + result.getMutationAA();
      if (!entries.containsKey(key)) {
        entries.put(key, new ArrayList<MutatedResult>());
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
            for (MutatedResult r : entries.get(key)) {
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
  private static void populateCodons(ArrayList<MutatedResult> results, Sequence dna) {
    
    for (MutatedResult result : results) {
      
      String entry = dna.getMatchingEntry(result.getProteinName());
      int shift = result.getOffset();
      int offset = result.getShift();
      if (!result.isReversed()) {
        if (!result.isInsertion()) {
          result.setCodon(entry.substring((int)(result.getMutationPosition()*3+offset-3), (int)(offset+result.getMutationPosition()*3)).toUpperCase());
          if (Codon.translate(result.getCodon())!=result.getOriginalAA()) {
            System.err.printf("Incorrect forward codons %c\t%s\t%s\t%d\n", result.getOriginalAA(), result.getCodon(), result.getProteinName(), offset);
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
          long position = end - shift - 3*result.getMutationPosition();
          StringBuffer sb = new StringBuffer(entry.substring((int)position, (int)(position+3)).toUpperCase()).reverse();
          StringBuffer comp = new StringBuffer();
          for (int i=0; i<sb.length(); i++) {
            comp.append(Codon.complement(sb.charAt(i)));
          }
          result.setCodon(comp.toString());
          if (Codon.translate(result.getCodon())!=result.getOriginalAA()) {
            System.err.printf("Incorrect reverse codons %c\t%s\t%d\t%d\t%s\n", result.getOriginalAA(), result.getCodon(), offset, shift, result.getProteinName());
          }
        }
      }
    }
  }
  
  
  /**
   * Keep only those result items that lie within a error tolerance from 
   * experimental to theoretical mass
   * @param results
   * @param cutOff
   * @return
   */
  private static ArrayList<MutatedResult> keepByDelta(ArrayList<MutatedResult> results, float cutOff) {
    ArrayList<MutatedResult> retItems = new ArrayList<MutatedResult>();
    for (MutatedResult result : results) {
      if (Math.abs(result.getMassError()) <= cutOff) retItems.add(result);
    }
    return retItems;
  }
  
  
  
  /**
   * Filter out all those items that do not have the mutation at the beginning
   * @param results the list of results
   * @return the filtered list of results
   */
  private static ArrayList<MutatedResult> keepFirstPosition(ArrayList<MutatedResult> results) {
    ArrayList<MutatedResult> retItems = new ArrayList<MutatedResult>();
    for (MutatedResult result : results) {
      if (result.getMutationPosition()==result.getStartPosition()) retItems.add(result);
    }
    return retItems;
  }
  
  
  private static ArrayList<MutatedResult> probFilter(ArrayList<MutatedResult> results, float prob) {
    ArrayList<MutatedResult> retItems = new ArrayList<MutatedResult>();
    for (MutatedResult result : results) {
      if (result.getProb() <= prob) retItems.add(result);
    }
    return retItems;
  }
  
  
  
  public static void processMultipleMutatedResultsAvar() {
    
    String userHome = System.getProperty("user.home");

    String[] resultFiles =  {userHome + "/Data/Spectra/Avar/ORG013_LTQ_0/allResultsMutated.txt"};
                           //  userHome + "/Data/Spectra/Avar/ORG013_LTQ_1/allResultsMutated.txt"};
    String dnaDbFile = userHome + "/Data/Databases/Avar/gen/Avar.fasta";
    String outDir = userHome + "/Data/Spectra/Avar/allResultsMutationTable";

    Sequence dnaDb = new FastaSequence(dnaDbFile);
    
    ArrayList<MutatedResult> results = new ArrayList<MutatedResult>();
    
    try {
      for (String resultFile : resultFiles) {
        BufferedReader in = new BufferedReader(new FileReader(resultFile));
        
        String line;
        while ((line = in.readLine())!=null) {
          MutatedResult current = new MutatedResult(line);
          results.add(current);
        }
      }
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    int initialSize = results.size();
    ArrayList<MutatedResult> filteredResults = results;
    filteredResults = probFilter(filteredResults, 1e-11f);
    filteredResults = keepFirstPosition(filteredResults);
    filteredResults = keepByDelta(filteredResults, 0.05f); 
    //filteredResults = eliminateRedundant(filteredResults);
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
    
    ArrayList<MutatedResult> results = new ArrayList<MutatedResult>();
    
    try {
      for (String resultFile : resultFiles) {
        BufferedReader in = new BufferedReader(new FileReader(resultFile));
        
        String line;
        while ((line = in.readLine())!=null) {
          if (line.split("\\s").length<14) continue;
          MutatedResult current = new MutatedResult(line);
          results.add(current);
        }
      }
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    int initialSize = results.size();
    ArrayList<MutatedResult> filteredResults = results;
    filteredResults = probFilter(filteredResults, 1e-13f);
    filteredResults = keepFirstPosition(filteredResults);
    filteredResults = keepByDelta(filteredResults, 0.02f); 
    //filteredResults = eliminateRedundant(filteredResults);
    populateCodons(filteredResults, dnaDb);
    makeMutationTable(outDir, filteredResults);
    System.out.printf("Previous count: %d. Filtered count: %d\n", initialSize, filteredResults.size());    
  }
  
  

  public static void processMultipleMutatedResultsAsp() {
    
    String userHome = System.getProperty("user.home");

    String[] resultFiles = {userHome + "/Data/Spectra/Asp/ORG014_LTQ_FT_0/results6Mutated.txt"};
                            //userHome + "/Data/Spectra/Asp/ORG014_LTQ_3/allResultsMutated.txt",
                            //userHome + "/Data/Spectra/Asp/ORG014_LTQ_2/allResultsMutated.txt",
                            //userHome + "/Data/Spectra/Asp/ORG014_LTQ_1/allResultsMutated.txt",
                            //userHome + "/Data/Spectra/Asp/ORG014_LTQ_0/allResultsMutated.txt"};
    String dnaDbFile = userHome + "/Data/Databases/Asp/gen/Asp.fasta";
    String outDir = userHome + "/Data/Spectra/Asp/mutationTable";

    Sequence dnaDb = new FastaSequence(dnaDbFile);
    
    ArrayList<MutatedResult> results = new ArrayList<MutatedResult>();
    
    try {
      for (String resultFile : resultFiles) {
        BufferedReader in = new BufferedReader(new FileReader(resultFile));
        
        String line;
        while ((line = in.readLine())!=null) {
          if (line.startsWith("#")) continue;
          MutatedResult current = new MutatedResult(line);
          results.add(current);
        }
      }
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    
    int initialSize = results.size();
    
    ArrayList<MutatedResult> filteredResults = results;
    filteredResults = probFilter(filteredResults, 1.50e-15f);
    filteredResults = keepFirstPosition(filteredResults);
    filteredResults = keepByDelta(filteredResults, 0.02f); 
    //filteredResults = eliminateRedundant(filteredResults);
    populateCodons(filteredResults, dnaDb);
    makeMutationTable(outDir, filteredResults);
    System.out.printf("Previous count: %d. Filtered count: %d\n", initialSize, filteredResults.size());    
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
    
    ArrayList<MutatedResult> results = new ArrayList<MutatedResult>();
    
    try {
      for (String resultFile : resultFiles) {
        BufferedReader in = new BufferedReader(new FileReader(resultFile));
        
        String line;
        while ((line = in.readLine())!=null) {
          MutatedResult current = new MutatedResult(line);
          results.add(current);
        }
      }
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    int initialSize = results.size();
    ArrayList<MutatedResult> filteredResults = results;
    filteredResults = probFilter(filteredResults, 1e-13f);
    filteredResults = keepFirstPosition(filteredResults);
    filteredResults = keepByDelta(filteredResults, 0.02f); 
    //filteredResults = eliminateRedundant(filteredResults);
    populateCodons(filteredResults, dnaDb);
    makeMutationTable(outDir, filteredResults);
    System.out.printf("Previous count: %d. Filtered count: %d\n", initialSize, filteredResults.size());    
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
    
    ArrayList<MutatedResult> results = new ArrayList<MutatedResult>();
    
    try {
      for (String resultFile : resultFiles) {
        BufferedReader in = new BufferedReader(new FileReader(resultFile));
        
        String line;
        while ((line = in.readLine())!=null) {
          if (line.split("\t").length < 12 ) continue;
          MutatedResult current = new MutatedResult(line);
          results.add(current);
        }
      }
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    int initialSize = results.size();
    ArrayList<MutatedResult> filteredResults = results;
    filteredResults = probFilter(filteredResults, 1e-11f);
    //filteredResults = keepFirstPosition(filteredResults);
    filteredResults = keepByDelta(filteredResults, 0.02f); 
    //filteredResults = eliminateRedundant(filteredResults);
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
    
    ArrayList<MutatedResult> results = new ArrayList<MutatedResult>();
    
    try {
      for (String resultFile : resultFiles) {
        BufferedReader in = new BufferedReader(new FileReader(resultFile));
        
        String line;
        while ((line = in.readLine())!=null) {
          if (line.split("\t").length < 12 ) continue;
          MutatedResult current = new MutatedResult(line);
          results.add(current);
        }
      }
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    int initialSize = results.size();
    ArrayList<MutatedResult> filteredResults = results;
    filteredResults = probFilter(filteredResults, 1e-13f);
    //filteredResults = keepFirstPosition(filteredResults);
    filteredResults = keepByDelta(filteredResults, 0.02f); 
    //filteredResults = eliminateRedundant(filteredResults);
    populateCodons(filteredResults, dnaDb);
    makeMutationTable(outDir, filteredResults);
    System.out.printf("Previous count: %d. Filtered count: %d\n", initialSize, filteredResults.size());    
  }


  public static void main(String[] args) {
    //processMutatedResults();
    
    processMultipleMutatedResultsAsp();
    //processMultipleMutatedResultsEcoli();
    //processMultipleMutatedResultsCsp();
    //processMultipleMutatedResultsAvar();
    //processMultipleMutatedResultsScerv();
    //processMultipleMutatedResultsScervMito();
    //processMultipleModdedResultsSone();
  }
}
