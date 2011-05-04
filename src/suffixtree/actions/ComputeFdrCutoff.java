package suffixtree.actions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;


/**
 * Given a forward and reverse search results compute the false discovery 
 * rate at a given cut off
 * @author jung
 *
 */
public class ComputeFdrCutoff {
  
  /**
   * Read the scores from the result files and store them into a sorted set.
   * Do no count duplicates
   * @param results the path to the result file
   * @return the set of sorted results
   */
  private static Float[] readScores(String results) {
    ArrayList<Float> scores = new ArrayList<Float>();
    try {
      BufferedReader in = new BufferedReader(new FileReader(results));
      String line, prevFile = "";
      int prevScanNum = -1;
      while ((line=in.readLine())!=null) {
        if (line.startsWith("#")) continue;
        
        String[] tokens = line.split("\t");
        int scanNum = Integer.parseInt(tokens[1]);
        String file = tokens[0];
        if (scanNum!=prevScanNum || !file.equals(prevFile)) {
          float prob = Float.parseFloat(tokens[6]);
          scores.add(prob);
        }
        prevScanNum = scanNum; prevFile = file;
      }
      
    }
    catch(IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    Collections.sort(scores);
    return scores.toArray(new Float[scores.size()]);
  }
  
  
  /**
   * Main method to compute the fdr given 2 result files
   * @param fResults the forward search result file
   * @param bResults the reverse search result file
   * @param percCutoff the fdr cutoff for the fdr
   */
  public static void computeFdr(String fResults, String bResults, float percCutoff) {
    
    Float[] fScores = readScores(fResults);
    Float[] rScores = readScores(bResults);
    System.out.printf("Total forward scores %d. Total reverse scores %d\n", fScores.length, rScores.length);
    
    float currentScore = Math.max(fScores[0], rScores[0]);
    int cumPos = 0, cumNeg = 0;
    int fIndex = 0, rIndex = 0;
    while (true) {
      while(fScores[fIndex] <= currentScore) {
        fIndex++;
        cumPos++;
      }
      
      while(rScores[rIndex] <= currentScore) {
        rIndex++;
        cumNeg++;
      }
      
      float fdr = (cumNeg * 100.0f) / (cumPos+cumNeg);
      System.out.printf("Good Spectra %d\tScore %.2e\tFdr %.2f\n", fIndex, currentScore, fdr);
      
      if (fdr > percCutoff) break;
      currentScore = Math.max(fScores[fIndex], rScores[rIndex]);
    }
  }

  
  
  public static void main(String[] args) {
    String fFile, rFile;
    float cutOff;
    
    if (args.length==0) {
      String userHome = System.getProperty("user.home");
      fFile = userHome + "/Data/Spectra/Sone/LTQFT0/bench.txt";
      rFile = userHome + "/Data/Spectra/Sone/LTQFT0/benchR.txt";
      
      fFile = userHome + "/Data/Spectra/Asp/ORG014_LTQ_FT_0/results6.txt";
      rFile = userHome + "/Data/Spectra/Asp/ORG014_LTQ_FT_0/results6R.txt";
      
      fFile = userHome + "/Data/Spectra/Asp/ORG014_LTQ_FT_0/results6Mutated.txt";
      rFile = userHome + "/Data/Spectra/Asp/ORG014_LTQ_FT_0/results6MutatedR.txt";
      
      //fFile = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/lys/output6.txt";
      //rFile = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/lys/output6R.txt";
      
      fFile = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/tryp/output6.txt";
      rFile = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/tryp/output6R.txt";
      
      fFile = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/tryp/output6Modded.txt";
      rFile = userHome + "/Data/Spectra/Hsapiens/Heck/mzXML/tryp/output6ModdedR.txt";
      
      cutOff = 1;
      computeFdr(fFile, rFile, cutOff);
    }
    else if (args.length==3) {
      fFile = args[0];
      rFile = args[1];
      cutOff = Float.parseFloat(args[2]);
      computeFdr(fFile, rFile, cutOff);
    }
    else {
      System.out.println("java ComputeFdrCutoff fowardResults reverseResults percCutOff");
    }
  }
  
  
  
}
