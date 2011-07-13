package scripts;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

import msutil.Composition;
import msutil.Peptide;
import msutil.Spectrum;

import parser.MzXMLSpectraMap;

public class SelectSpectra {

  public static void select() {
    String userHome = System.getProperty("user.home");
    
    String idFile = userHome + "/Desktop/PAe000353_mzXML_200903080810/results/interact-combined.pep.PAidentlist";
    String specDir = userHome + "/Desktop/PAe000353_mzXML_200903080810";
    String outFile = userHome + "/Desktop/MmusHeartMito.ms2";
    String ids = userHome + "/Desktop/MmusHeartMitoIds.fasta";
    
    //String idFile = userHome + "/Desktop/PAe000359_mzXML_200903071001/results/interact-combined.pep.PAidentlist";
    //String specDir = userHome + "/Desktop/PAe000359_mzXML_200903071001";
    //String outFile = userHome + "/Desktop/MmusBrainMito.ms2";
    //String ids = userHome + "/Desktop/MmusBrainMitoIds.fasta";
    
    String ipiDb = userHome + "/Desktop/ipi.MOUSE.v3.74.fasta";
    HashMap<String,String> id2desc = new HashMap<String,String>();
    
    // parse the ipi file
    try {
      BufferedReader in = new BufferedReader(new FileReader(ipiDb));
      String line = null;
      while ((line = in.readLine())!=null) {
        if (line.startsWith(">")) {
          String ipiId = line.split("IPI", 2)[1].substring(1).split(":", 2)[0].split("\\.")[0];
          String anno = line.split("Gene_Symbol", 2)[1].substring(1);
          id2desc.put(ipiId, anno);
          //System.out.println(ipiId + "\t" + anno);
        }
      }
      
    }
    catch (IOException ioe) {
      System.err.println(ioe);
    }
    
    int totalCount = 0;
    int mitoPro = 0;
    
    try {
      BufferedReader in = new BufferedReader(new FileReader(idFile));
      PrintWriter fout = new PrintWriter(outFile);
      PrintWriter idOut = new PrintWriter(ids);
      
      in.readLine(); // skip first line
      
      String currentFile = "";
      MzXMLSpectraMap sm = null;
      
      String line = in.readLine();
      while (line!=null) {
        String[] tokens = line.split("\t");
        String[] infoTokens = tokens[1].split("\\.");
        
        String filename = infoTokens[0]; 
        int scanNum = Integer.parseInt(infoTokens[1]);
        int charge = Integer.parseInt(infoTokens[3]);
        String ipiId = tokens[10];
        
        String anno = id2desc.get(ipiId);
        if (anno==null) {
          //System.out.println(ipiId + " not found");
        }
        else {
          //System.out.println(anno + " match!");
          if (anno.matches(".*itochondria.*")) {
            System.out.println(anno);
            mitoPro++;
          }
        }
        
        //float massDiff = Float.parseFloat(tokens[9]);
        
        if (!filename.equals(currentFile)) {
          sm = new MzXMLSpectraMap(specDir + "/" + filename + ".mzXML");
          currentFile = filename;
        }
        
        Spectrum s = sm.getSpectrumBySpecIndex(scanNum);
        s.setCharge(charge);
        
        //System.out.printf("%.2f %.2f\n", s.getParentMass(), new Peptide(tokens[3]).getMass());
        
        // exact mass
        float pm = s.getParentMass();
        double em = new Peptide(tokens[3]).getMass() + Composition.H2O;// + Composition.H;
        float massCorrection = (float)em-pm;
        s.setPrecursor(s.getPrecursorPeak().duplicate(massCorrection));
        
        fout.printf(":%d.%d.0\n", totalCount, totalCount);
        fout.println(s.toDta());
        
        idOut.printf(">%d\n%s\n", totalCount, tokens[4]+tokens[3]+tokens[6]);

        totalCount++;
        
        if (totalCount > 1000) break;
        
        line = in.readLine();
      }
      
      fout.close();
      in.close();
      idOut.close();
      
    }
    catch (IOException ioe) {
      System.err.println(ioe);
    }
    
    System.out.println("Total spec " + totalCount + ". Total mito proteins " + mitoPro);
  }
  
  
  public static void main(String[] args) {
    select();  
  }
}
