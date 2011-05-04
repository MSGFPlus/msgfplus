package suffixtree.actions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Iterator;

import msgap.results.GappedPeptideResults;
import msgap.results.MSGDResultFileParser;

/**
 * This class allows one process the grc file so that only the unmatched 
 * spectra's gapped peptide are retained
 * @author jung
 *
 */
public class TrimGRCFile {

  
  public static void retain(String resultFile, float cutoff, String grcIn, String grcOut) {
    
    HashSet<String> goodMatches = new HashSet<String>();
    PrintWriter out = null;
    
    // register all passing matches 
    try {
      
      BufferedReader in = new BufferedReader(new FileReader(resultFile));
      String line;
      while ((line=in.readLine())!=null) {
        if (line.startsWith("#")) continue;
        String[] tokens = line.split("\t");
        float prob = Float.parseFloat(tokens[6]);
        if (prob <= cutoff) {
          String key = tokens[1] + "%%" + tokens[0];    // scan%file is the key
          goodMatches.add(key);
        }
      }
      
      // open the output grc file
      out = new PrintWriter(grcOut);
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    // calculate how many results to parse at a time based on the total free memory
    int queryCount = (int)(Runtime.getRuntime().maxMemory() * 0.00008);
    Iterator<GappedPeptideResults> it = new MSGDResultFileParser(grcIn, queryCount).iterator();
    int totalGpr = 0;
    int leftGpr = 0;
    while (it.hasNext()) {
      HashSet<Integer> goodIds = new HashSet<Integer>(); 
      GappedPeptideResults gpr = it.next();
      totalGpr += gpr.getSpectrumCount();
      
      for (int specId : gpr.getSpecIds()) {
        String key = gpr.getScanNumber(specId) + "%%" + gpr.getFileName(specId);
        if (!goodMatches.contains(key)) {
          goodIds.add(specId);
        }
      }
      
      GappedPeptideResults leftovers = gpr.retain(goodIds);
      leftGpr += leftovers.getSpectrumCount();
      leftovers.toFile(out);
      
    }
    
    System.out.printf("Total Spectra %d\tSpectra with no good matches %d\tProb cutoff %.2e\n", totalGpr, leftGpr, cutoff);
    
  }
  
  
  public static void main(String[] args) {
    String rFile, grcIn, grcOut;
    float cutOff;
    
    if (args.length==0) {
      String userHome = System.getProperty("user.home");
      rFile = userHome + "/Data/Spectra/Sone/LTQFT0/bench.txt";
      cutOff = 3.5e-12f;
      grcIn = userHome + "/Data/Spectra/Sone/LTQFT0/bench.grc";
      grcOut = userHome + "/Data/Spectra/Sone/LTQFT0/benchMutated.grc";
      
      rFile = userHome + "/Data/Spectra/Asp/ORG014_LTQ_FT_0/results6.txt";
      cutOff = 1.70e-10f;
      grcIn = userHome + "/Data/Spectra/Asp/ORG014_LTQ_FT_0/results6.grc";
      grcOut = userHome + "/Data/Spectra/Asp/ORG014_LTQ_FT_0/results6Mutated.grc";
      
      /*
      rFile = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/lys/output6.txt";
      cutOff = 6.2e-12f;
      grcIn = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/lys/output6.grc";
      grcOut = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/lys/output6Mutated.grc";
      */
      
      rFile = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/tryp/output6.txt";
      cutOff = 4.00e-11f;
      grcIn = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/tryp/output6.grc";
      grcOut = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/tryp/output6Mutated.grc";
      
      retain(rFile, cutOff, grcIn, grcOut);
    }
    else if (args.length==4) {
      rFile = args[0];
      cutOff = Float.parseFloat(args[1]);
      retain(rFile, cutOff, args[2], args[3]);
    }
    else {
      System.out.println("java TrimGRCFile resultFile probCutOff grcIn grcOut");
    }
  }
  
}
