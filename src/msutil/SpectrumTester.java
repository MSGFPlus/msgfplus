/**
 * Utilities.
 */
package msutil;

import java.io.*;
import java.util.*;


/**
 * This class tests the Spectrum and its dependent classes.
 * @author jung
 *
 */
public class SpectrumTester {

  /**
   * Main.
   * @param args main parameters.
   */
  public static void main(String[] args) {
    long time = System.currentTimeMillis();
    String mgfFile = "/home/sangtaekim/Research/Data/PNNL/IPYS_TD_Scere010_Orbitrap_001a.mgf";
    parseMgf(mgfFile);
    System.out.println("Time: " + (System.currentTimeMillis() - time));
  }

  // temp method to parse an MGF file from a path
  public static ArrayList<Spectrum> parseMgf(String path) {
    
    BufferedReader fileReader = null;
    ArrayList<Spectrum> results = new ArrayList<Spectrum>();
    try {
      fileReader = new BufferedReader(new FileReader(path));
      boolean parseFlag = false;
      String line = null;
      Spectrum tempSpec = null;        // the spectrum being parsed       
      int specCount = 0;
      while((line = fileReader.readLine()) != null) {
        
        // try to parse a peak
        if(parseFlag && Character.isDigit(line.charAt(0))) {
          // parse as a peak
          String[] tokens = line.split("\\s");
          float intensity = Float.parseFloat(tokens[1]);
          if(intensity>Constants.EPSILON) {
            // only include non-zero peaks
            tempSpec.add(new Peak(Float.parseFloat(tokens[0]), intensity, 1));
          }
        }
        // try to parse the keywords
        else {
          String[] tokens = line.split("=");
          String keyword = tokens[0];
          if(keyword.equals("BEGIN IONS") && !parseFlag) {
            tempSpec = new Spectrum();
            parseFlag = true;
          }
          else if(parseFlag) {
            if(keyword.equals("END IONS")) {
              // create the spectrum object-*
              results.add(tempSpec);
              parseFlag = false;
              specCount += 1;
//              System.out.println(specCount + " " + tempSpec.size());
            }
            else if(keyword.equals("PEPMASS")) {
              float parent = Float.parseFloat(tokens[1]);
              if(tempSpec.getPrecursorPeak() != null) {
                tempSpec.getPrecursorPeak().setMz(parent);
              }
              else {
                tempSpec.setPrecursor(new Peak(parent, 1, 1));
              }
            }
            else if(keyword.equals("CHARGE")) {
              // the charge string is a little funky when it has a trailing charge 
              String chargeStr = tokens[1];
              int lastIndex = chargeStr.length()-1;
              char lastChar = chargeStr.charAt(lastIndex);
              if(lastChar == '-') {
                chargeStr = lastChar + chargeStr.substring(0, lastIndex);
              }
              else if(lastChar == '+') {
                chargeStr = chargeStr.substring(0, lastIndex);
              }
              if(chargeStr.charAt(0)=='+') {
                chargeStr = chargeStr.substring(1); 
              }
              int charge = java.lang.Integer.parseInt(chargeStr);
              if(tempSpec.getPrecursorPeak() != null) {
                tempSpec.getPrecursorPeak().setCharge(charge);
              }
              else {
                tempSpec.setPrecursor(new Peak(0.0f, 1, charge));
              }
            }
            else if(keyword.equals("TITLE")) {
              tempSpec.setTitle(tokens[1]);
            }
          }
        }
      }
    }
    catch(IOException e) { 
      System.err.println(e.getMessage());
      return null;
    }
    
    
    
  return results;  
  }
}
