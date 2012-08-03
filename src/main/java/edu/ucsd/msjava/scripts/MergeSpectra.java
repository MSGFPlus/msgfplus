package edu.ucsd.msjava.scripts;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.MS2SpectrumParser;
import edu.ucsd.msjava.parser.MgfSpectrumParser;
import edu.ucsd.msjava.parser.MzXMLSpectraIterator;


/**
 * Merge a certain count of spectra into a single file.
 * @author jung
 *
 */
public class MergeSpectra {

  public static void mergeMs2IntoMs2(String inputDir, String outFile, int count) {
    
    PrintWriter fout = null;
    try {
      fout = new PrintWriter(outFile);
    }
    catch(IOException ioe) {
      System.err.println(ioe);
      System.exit(-9);
    }
    
    int totalCount = 0;
    File dir = new File(inputDir);
    for (String fileName : dir.list()) {
      try {
        SpectraIterator sIt = new SpectraIterator(inputDir + "/" + fileName, new MS2SpectrumParser());
        while (sIt.hasNext()) {
          fout.printf(":%d.%d.0\n", totalCount, totalCount);
          fout.println(sIt.next().toDta());
          totalCount++;
          if (totalCount>=count) {
            fout.close();
            return;
          }
        }
      }
      catch (IOException ioe) {
        System.err.println(ioe);
        System.exit(-1);
      }
    }
  }

  
 public static void mergeMgfIntoMs2(String inputDir, String outFile, int count) {
    
    PrintWriter fout = null;
    try {
      fout = new PrintWriter(outFile);
    }
    catch(IOException ioe) {
      System.err.println(ioe);
      System.exit(-9);
    }
    
    int totalCount = 0;
    File dir = new File(inputDir);
    for (String fileName : dir.list()) {
      try {
        if (!fileName.endsWith(".mgf")) continue;
        SpectraIterator sIt = new SpectraIterator(inputDir + "/" + fileName, new MgfSpectrumParser());
        while (sIt.hasNext()) {
          
          Spectrum s = sIt.next();
          if (s.getParentMass() < 10) continue;
        
          fout.printf(":%d.%d.0\n", totalCount, totalCount);
          fout.println(s.toDta());
          totalCount++;
          if (totalCount>=count) {
            fout.close();
            return;
          }
        }
      }
      catch (IOException ioe) {
        System.err.println(ioe);
        System.exit(-1);
      }
    }
  }
  
  

 
 public static void mergeMzXMLIntoMs2(String inputDir, String outFile, int count) {
    
    PrintWriter fout = null;
    try {
      fout = new PrintWriter(outFile);
    }
    catch(IOException ioe) {
      System.err.println(ioe);
      System.exit(-9);
    }
    
    int totalCount = 0;
    File dir = new File(inputDir);
    for (String fileName : dir.list()) {
      MzXMLSpectraIterator sIt = new MzXMLSpectraIterator(inputDir + "/" + fileName); 
      while (sIt.hasNext()) {
        Spectrum s = sIt.next();
        
        // get enough peaks
        if (s.size() < 30) continue;
        
        s.setCharge(2);
        fout.printf(":%d.%d.0\n", totalCount, totalCount);
        fout.println(s.toDta());
        totalCount++;
        if (totalCount>=count) {
          fout.close();
          return;
        }
      }
    }
  }
  
 
  
  public static void main(String[] args) {
    String userHome = System.getProperty("user.home");
    
    int totalSpec = 10000000;
    //String inputDir = userHome + "/Data/Spectra/Sone/LTQFT0";
    //String inputDir = userHome + "/Data/Spectra/SCerv/select";
    //String inputDir = userHome + "/Desktop/PAe000359_mzXML_200903071001";
    String inputDir = "/home/jung/Data/Spectra/Maize/SQS02/Maize-Pericarp-Aleurone";
    
    //String outputFile = userHome + String.format("/Desktop/Sone%d.ms2", totalSpec);
    String outputFile = userHome + String.format("/Desktop/Maize-Pericarp-Aleurone3.ms2");
    
    //mergeMs2IntoMs2(inputDir, outputFile, totalSpec);
    mergeMgfIntoMs2(inputDir, outputFile, totalSpec);
    //mergeMzXMLIntoMs2(inputDir, outputFile, totalSpec);
  }
}
