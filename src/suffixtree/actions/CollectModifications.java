package suffixtree.actions;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeMap;

import suffixtree.results.ModdedResult;
import cyclic.Cluster1D;
import cyclic.Point1D;

public class CollectModifications {

  
  
  public static void collectOffsetDist(String[] dirs, float probCutoff, String outdir) {
    
    
    ArrayList<ModdedResult> results = new ArrayList<ModdedResult>();
    HashSet<String> keys = new HashSet<String>();
    try {
      for (String resultDir : dirs) {
        BufferedReader in = new BufferedReader(new FileReader(resultDir));
        String line;
        while ((line=in.readLine())!=null) {
          if (line.startsWith("#")) continue;
          
          ModdedResult r = new ModdedResult(line);
          String key = r.getFilepath() + "$$" + r.getScanNumber();
          if (!keys.contains(key)) {
            if (r.getProb()<=probCutoff) results.add(r);           
            keys.add(key);
          }
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
    for (ModdedResult r : results) {
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
      
      // create the directory
      if (!new File(outdir + "/Mods").exists()) {
        new File(outdir + "/Mods").mkdirs();
      }
      
      out.println("<table border='1'><tr><th>Offset</th><th>Spectra</th><th>Peptides</th></tr>");
      for (Cluster1D c : clusters) {
        HashSet<String> sites = new HashSet<String>();
        // print out to a separate file
        String filename = String.format("Mods/file_%.3f.html", c.getCenter());
        String filepath = String.format("%s/%s", outdir, filename);
        PrintWriter pw = new PrintWriter(filepath);
        pw.println("<table border='1'><tr><th>Peptide</th><th>Delta</th><th>Charge</th><th>Filename</th><th>Scan</th><th>Prob</th></tr>");
        for (Point1D p : c.getPoints()) {
          ModdedResult r = resultArray.get(p.getIndex());
          sites.add(r.getModificationPositionKey());
          pw.printf("<tr><td>%s</td><td>%.3f</td><td>%d</td><td>%s</td><td>%d</td><td>%.3e</td></tr>", r.getPeptide(), r.getDelta(), r.getCharge(), r.getFilename(), r.getScanNumber(), r.getProb());
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
  
  
  public static void processMultipleModdedResultsHuman() {
    String userHome = System.getProperty("user.home");
    
    String[] resultFiles = {userHome+"/Data/Spectra/Hsapiens/Heck/mzXML/tryp/output6Modded.txt"};
    
  
    // print out the offset counts given the cutoff
    String offsetDir = userHome+"/Data/Spectra/Hsapiens/Heck/mzXML/tryp";
    
    collectOffsetDist(resultFiles, 1.30e-15f, offsetDir);
  }
  
  
  public static void main(String[] args) {
    processMultipleModdedResultsHuman();
  }
  
  
}
