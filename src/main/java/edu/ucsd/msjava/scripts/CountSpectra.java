package edu.ucsd.msjava.scripts;

import java.io.File;
import java.io.IOException;

import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.parser.MS2SpectrumParser;
import edu.ucsd.msjava.parser.MzXMLSpectraIterator;

/**
 * Count the spectra in a directory
 * @author jung
 *
 */
public class CountSpectra {

  public static int countMS2(String path) {
    File fObj = new File(path);
    if (fObj.isDirectory()) {
      int count = 0;
      String[] files = fObj.list(); 
      for (String file : files) {
        if (file.endsWith(".ms2")) {
          count += countMS2file(path + "/" + file);
        }
      }
      return count;
    }
    else {
      return countMS2file(path);
    }
  }
  
  private static int countMS2file(String file) {
    try {
      int count = 0;
      SpectraIterator it = new SpectraIterator(file, new MS2SpectrumParser());
      while (it.hasNext()) {
        it.next();
        count += 1;
      }
      return count;
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    return 0;
  }
  
  public static int countMZXML(String path) {
    File fObj = new File(path);
    if (fObj.isDirectory()) {
      int count = 0;
      String[] files = fObj.list(); 
      for (String file : files) {
        if (file.endsWith(".mzXML")) {
          count += countMZXMLfile(path + "/" + file);
        }
      }
      return count;
    }
    else {
      return countMS2file(path);
    }
  }
  
  private static int countMZXMLfile(String file) {
    int count = 0;
    MzXMLSpectraIterator it = new MzXMLSpectraIterator(file);
    while (it.hasNext()) {
      it.next();
      count += 1;
    }
    return count;
  }
  
  
  public static void main(String[] args) {
    String userHome = System.getProperty("user.home");
    
    //String dir = String.format("%s/Data/Spectra/Sone/LTQFT0", userHome);
    //int count = countMS2(dir);
    
    String dir = String.format("%s/Data/Spectra/Hsapiens/Heck/mzXML/tryp", userHome);
    int count = countMZXML(dir);
    
    System.out.println("Number of spectra: " + count);
  }
  
  
}
