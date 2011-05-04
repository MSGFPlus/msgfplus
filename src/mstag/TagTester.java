package mstag;

import msutil.*;
import java.io.*;

import parser.MgfSpectrumParser;

/**
 * Testing code for high accuracy tags
 * @author jung
 *
 */
public class TagTester {

  /**
   * @param args
   */
  public static void main(String[] args) throws FileNotFoundException {
    String filePath = "/Users/jung/Research/USTags/centroidedMGF/IPYS_TD_Scere010_Orbitrap_001c.mgf";
    SpectraIterator si = new SpectraIterator(filePath, new MgfSpectrumParser());

    Reshape windowFilter = new WindowFilter(5, 50);
    int specCount = 0;
    while(si.hasNext()) {
      Spectrum spec = si.next();
      System.out.print("Raw peak count: "+spec.size()+'\t');
      spec = windowFilter.apply(spec);
      System.out.println("Filtered peak count: "+spec.size());
      if(spec.getAnnotationStr()!=null) {
        System.out.println(spec.getAnnotationStr());
        specCount++;
      }
      
    }
    System.out.println("Size of container is " + specCount);

  }

}
