package suffixtree.misc;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

public class MatchingStats {

  private static class Stat {
    private float prob, offset;
    
    private Stat(float prob, float offset) {
      this.prob = prob;
      this.offset = offset;
    }
  }
  
  // holds the scores gotten per each match of the spectra
  private HashMap<Integer,ArrayList<Stat>> stats;
  private int statCount;
  
  public MatchingStats() {
    this.stats = new HashMap<Integer,ArrayList<Stat>>();
    this.statCount = 0;
  }
  
  
  public void addItem(int specId, float prob, float offset) {
    if (!this.stats.containsKey(specId)) {
      this.stats.put(specId, new ArrayList<Stat>());
    }
    this.stats.get(specId).add(new Stat(prob, offset));
    statCount++;
  }

  
  public void printOutScoreDist(PrintWriter out) {
    
    out.println("\n---> Score distribution");
    
    // the histogram using the log of the score
    TreeMap<Integer,ArrayList<Float>> scoreDist = new TreeMap<Integer,ArrayList<Float>>();
    
    for (int specId : this.stats.keySet()) {
      for (Stat s : this.stats.get(specId)) {
        int key = (int)Math.log10(s.prob);
        if (!scoreDist.containsKey(key)) scoreDist.put(key, new ArrayList<Float>());
        scoreDist.get(key).add(s.prob);
      }
    }
    
    for (int key : scoreDist.keySet()) {
      out.printf("%d\t%d\n", key, scoreDist.get(key).size());
    }
  }
  
  
  public void printOutDualScoreDist(PrintWriter out) {
    
    out.println("\n---> Score distribution");
    
    // the histogram using the log of the score
    TreeMap<Integer,ArrayList<Float>> negScoreDist = new TreeMap<Integer,ArrayList<Float>>();
    TreeMap<Integer,ArrayList<Float>> posScoreDist = new TreeMap<Integer,ArrayList<Float>>();
    
    for (int specId : this.stats.keySet()) {
      for (Stat s : this.stats.get(specId)) {
        int key = (int)Math.log10(s.prob);
        if (s.offset >= 0) {
          if (!posScoreDist.containsKey(key)) posScoreDist.put(key, new ArrayList<Float>());
          posScoreDist.get(key).add(s.prob);
        }
        else {
          if (!negScoreDist.containsKey(key)) negScoreDist.put(key, new ArrayList<Float>());
          negScoreDist.get(key).add(s.prob);
        }
      }
    }

    out.printf("Positive offset score distribution");
    for (int key : posScoreDist.keySet()) {
      out.printf("%d\t%d\n", key, posScoreDist.get(key).size());
    }
    out.printf("Negative offset score distribution");
    for (int key : negScoreDist.keySet()) {
      out.printf("%d\t%d\n", key, negScoreDist.get(key).size());
    }
  }
  
  
  public void printMatchStats(PrintWriter out, float cutOff) {
    
    out.printf("\n---> Stats for %e prob cutoff\n", cutOff);
    
    // the histogram of offset distribution. Use 0.05Da
    float binSize = 0.01f;
    TreeMap<Integer,ArrayList<Float>> offsetDist = new TreeMap<Integer,ArrayList<Float>>();
    
    int specCount = 0, queryCount = 0;
    for (int specId : this.stats.keySet()) {
      boolean hasMatch = false;
      for (Stat s : this.stats.get(specId)) {
        if (s.prob <= cutOff) {
          hasMatch = true;
          queryCount++;
          
          int key = (int)((s.offset + binSize*0.5f) / binSize);
          if (!offsetDist.containsKey(key)) offsetDist.put(key, new ArrayList<Float>());
          offsetDist.get(key).add(s.offset);
        }
      }
      if (hasMatch) specCount++;
    }
    
    for (int key : offsetDist.keySet()) {
      out.printf("%.2f\t%d\n", key * binSize, offsetDist.get(key).size());
    }
    out.printf("%d spectra with %d matches.\n", specCount, queryCount);
  }
  
  
  public int getStatCount() { return statCount; }
  
  public void clear() { this.stats.clear(); }
  
}
