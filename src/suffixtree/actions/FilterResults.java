package suffixtree.actions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;


public class FilterResults {

  private final static float probCutOff = 1e-10f;
  
  
  public static boolean filterByProb(String[] tokens, float cutOff) {
    float p = Float.parseFloat(tokens[8]); 
    if (p <= cutOff && p != 0)    return true;
    return false;
  }
  
  public static boolean filterByCoors(String[] tokens, int lower, int upper) {
    int start = Integer.parseInt(tokens[9]);
    int end = Integer.parseInt(tokens[10]);
    if (lower < start && end < upper) return true;
    return false;
  }
  
  
  public static void processResultsCsp() {
    
    int coors = 685000;
    
    //coors = 89300;

    String userHome = System.getProperty("user.home");

    String[] resultFiles = {userHome + "/Data/Spectra/Csp/ORG033_LTQ_Orb_0/results6.txt",
                            userHome + "/Data/Spectra/Csp/ORG033_LTQ_0/results.txt"};
    String outFile = userHome + "/Data/Spectra/Csp/ORG033_LTQ_Orb_0/filtered.txt";
    
    int results = 0;
    
    try {
      PrintWriter out = new PrintWriter(outFile);
      
      for (String resultFile : resultFiles) {
        BufferedReader in = new BufferedReader(new FileReader(resultFile));
        
        String line;
        while ((line = in.readLine())!=null) {
          String[] tokens = line.split("\t");
          
          if (filterByProb(tokens, probCutOff) && filterByCoors(tokens, coors, coors+500)) {
            out.println(line);
            results++;
          }
          
        }
        
        in.close();
      }
      out.close();
    }
    catch (IOException ioe) {
      System.err.println(ioe);
      System.exit(-1);
    }
    
    System.out.printf("Matches %d\n", results);
    
  }
  
  
  public static void main(String[] args) {
    processResultsCsp();
  }
  
  
}
