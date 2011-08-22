package cyclic;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.ListIterator;
import java.util.Set;

import parser.MzXMLSpectraIterator;
import msutil.Spectrum;
import msutil.WindowFilter;


/**
 * Main class for the NRP-Tagging algorithm.
 * @author jung
 *
 */
public class Tag {
  
  public static void generateTriplets(Spectrum s, float tolerance) {
    
    // do double convolution to get the 3 masses
    Set<Cluster> clusters = Cluster.make2DClusters(s, Constants.MIN_DISTANCE, tolerance, Constants.MIN_CLUSTER_SIZE);
    
    // generate the sequence triplets
    
  }
  
  private static void debug() {
    String file;
    String home = System.getProperty("user.home");
    
    file = "/Data/Cyclic/spectra/TOF/tyrB1.mzXML";
    MzXMLSpectraIterator it = new MzXMLSpectraIterator(home + file, 1, 6);
   
    CircularSequence sequence = CircularSequence.tyrB1;
    float tolerance = 0.1f;
    
    WindowFilter filter = new WindowFilter(10, 50);
    
    Spectrum spec = null;
    while (it.hasNext()) {
      spec = it.next();
      spec = filter.apply(spec);
      spec.setPrecursorCharge(1);
    }
    
    // check how the original spectrum stacks up in correctness
    float[] peaks = new float[spec.size()];
    for (int i=0; i<peaks.length; i++) peaks[i] = spec.get(i).getMass();
    ArrayList<Float> matches = sequence.isSubset(peaks, 2*tolerance);
    System.out.println(spec);
    System.out.println(sequence);
    System.out.printf("Matches (%d): ", matches.size());
    for (float mass : matches) {
      System.out.printf("%.2f ", mass);
    }
    System.out.println();
    
    // do the 1D convolution and check how many correct masses are recovered
    Set<Cluster1D> guessedMasses = Cluster1D.make1DClusters(spec, Constants.MIN_DISTANCE, tolerance, Constants.MIN_CLUSTER_SIZE);
    //System.out.println("Number of guessed masses " + guessedMasses.size());
    
    ArrayList<Cluster1D> massList = new ArrayList<Cluster1D>();
    for (Cluster1D c : guessedMasses) {
      if (Constants.MIN_DISTANCE<c.getCenter() && c.getCenter()<Constants.MAX_DISTANCE) {
        massList.add(c);
      }
    }
    Collections.sort(massList, Cluster1D.weightComparator);
    if (massList.size() > Constants.MAX_CONV_MASSES) massList.subList(Constants.MAX_CONV_MASSES, massList.size()).clear();
    
    int clusterRank = 0;
    for (Cluster1D mass : massList) {
      clusterRank++;
      float[] masses = {mass.getCenter()};
      if (sequence.isSuperset(masses, tolerance)) {
        System.out.println(clusterRank + ": Correct " + mass.getCenter());
      }
      else {
        //System.out.println(clusterRank + ": Incorrect " + mass.getCenter());
      }
    }
    // convert the valid masses into an array
    float[] masses = new float[Constants.MAX_CONV_MASSES];
    for (int i=0; i<masses.length; i++) {
      masses[i] = massList.get(i).getCenter();
    }
    Arrays.sort(masses);
    
    // make the tags
    Set<Cluster> clusters = Cluster.make2DClusters(spec, Constants.MIN_DISTANCE, tolerance, Constants.MIN_CLUSTER_SIZE);
    System.out.println("Cluster count " + clusters.size());
    System.out.println();
    
    ArrayList<Cluster> clusterList = new ArrayList<Cluster>(clusters);
    Collections.sort(clusterList, new Cluster.ClusterWeightComparator());

    for (int i=0; i<30; i++) {
      generate(clusterList.get(i), spec, tolerance, sequence, masses);
      System.out.println(clusterList.get(i) + " ; " + sequence.isSuperset(clusterList.get(i).getCenters(), 2*tolerance));
      System.out.println();
    }

  }
  
  
  /**
   * Command line invoker.
   * @param args
   */
  public static void main(String[] args) {
    // testing
    debug();
  }

  
  
  /**
   * Generate the aligned spectrum.
   * @param c
   */
  private static void generate(Cluster c, Spectrum s, float tolerance, CircularSequence peptide, float[] masses) {
    ArrayList<Point> members = new ArrayList<Point>(c.getPoints());
    
    // calculate the offset for each point
    ArrayList<Float> offsets = new ArrayList<Float>(members.size());
    for (Point p : members) {
      float offset = 0.0f;
      for (int i=0; i<p.size()-1; i++) {
        offset += c.getCenters()[i] - p.getMassAt(i) + offset;
      }
      offsets.add((offset / p.size()) + s.get(p.getFirstIndex()).getMass());
    }
    
    // create a rotated version of the spectrum by shifting back by the offset value
    ArrayList<Point1D> pointList = new ArrayList<Point1D>();
    for (int i=0; i<members.size(); i++) {
      for (int j=0; j<s.size(); j++) {
        float shiftedMass = s.get(j).getMass()-offsets.get(i);
        if (shiftedMass < -2*tolerance) shiftedMass += s.getParentMass();
        if (shiftedMass <= s.getParentMass()-Constants.MIN_DISTANCE) {
          pointList.add(new Point1D(shiftedMass, s.get(j).getIntensity(), j));
        }
      }
    }
    
    // cheat the elimination algorithm by adding many points at the end
    for (int i=0; i<Constants.MIN_PEAK_COUNT; i++) pointList.add(new Point1D(s.getParentMass(), 1, 0));
    
    
    System.out.println("Total point count: " + pointList.size());
    
    Set<Cluster1D> positions = Cluster1D.cluster(pointList, 2*tolerance, Constants.MIN_PEAK_COUNT);
    System.out.println("Total cluster positions: " + positions.size());

    ArrayList<Cluster1D> clusterList = new ArrayList<Cluster1D>(positions);
    Collections.sort(clusterList);
    System.out.println("First point " + clusterList.get(0).getCenter());
    System.out.println("Last point " + clusterList.get(clusterList.size()-1).getCenter());
    
    float[] peaks = new float[positions.size()];
    int i = 0;
    for (Cluster1D cl : clusterList) {
      peaks[i++] = cl.getCenter();
    }
    
    Arrays.sort(peaks);
    ArrayList<Float> matches = peptide.isSubset(peaks, 2*tolerance);
    System.out.println("Best matching score " + peptide.getScore(clusterList, tolerance));
    System.out.printf("Matches (%d): ", matches.size());
    for (float mass : matches) {
      System.out.printf("%.2f ", mass);
    }
    System.out.println();
    
    // generate the best score
    //getBestScore(clusterList, Constants.MIN_DISTANCE, Constants.MAX_DISTANCE, masses, tolerance);
  } 
  
  
 /********** The code below is for sequencing, which doesn't really work **********/ 
  private static class Cell {
    private float bestScore;
    private float mass;
    private float score;
    private ArrayList<Cell> prevCells;
    private ArrayList<Float> scores;
    
    private Cell(float mass, float score) {
      this.mass = mass;
      this.score = score;
      this.prevCells = new ArrayList<Cell>();
      this.scores = new ArrayList<Float>();
      this.bestScore = 0;
    }
    
    private void addCell(Cell previous, float score) {
      this.prevCells.add(previous);
      this.scores.add(score+previous.bestScore);
      this.bestScore = Math.max(this.bestScore, score+previous.bestScore);
    }
  }

  
  private static void getBestScore(ArrayList<Cluster1D> points, float minDistance, float maxDistance, float[] masses, float tolerance) {
    
    ArrayList<Cell> table = new ArrayList<Cell>();
    table.add(new Cell(points.get(0).getCenter(), points.get(0).getWeight()));
    table.get(0).bestScore = points.get(0).getWeight();
    
    ListIterator<Cluster1D> outerIterator = points.listIterator(1);
    //System.out.println(points.get(0).getCenter() + " " + points.get(1).getCenter() + " " + outerIterator.next().getCenter());
    while (outerIterator.hasNext()) {
      Cluster1D point = outerIterator.next();
      
      Cell current = new Cell(point.getCenter(), point.getWeight());
      
      ListIterator<Cell> innerIterator = table.listIterator(table.size());
      while(innerIterator.hasPrevious()) {
        Cell previous = innerIterator.previous();
        
        if (FuzzyArray.search(masses, current.mass-previous.mass, tolerance)<0) continue; 
        if (current.mass - previous.mass > maxDistance) break;
        
        current.addCell(previous, point.getWeight());
      }
      
      table.add(current);
    }
    
    float bestScore = table.get(table.size()-1).bestScore;
    System.out.println("Best overall score: " + bestScore);
    
    ArrayList<ArrayList<Float>> results = generateSequences(table, bestScore*0.8f);
    System.out.println("Total Sequences generated " + results.size());
    
  }
  
  /**
   * Backtrack and generate all sequences with the given cutoff.
   * @param table the DP table with the scores
   * @param scoreCutoff the score cutoff threshold
   */
  private static ArrayList<ArrayList<Float>> generateSequences(ArrayList<Cell> table, float scoreCutoff) {
    ArrayList<ArrayList<Float>> results = new ArrayList<ArrayList<Float>>();
    Cell lastCell = table.get(table.size()-1);
    for (int i=0; i<lastCell.scores.size(); i++) {
      if (lastCell.scores.get(i)>=scoreCutoff) {
        //System.out.println("Score: " + lastCell.scores.get(i));
        ArrayList<Float> suffix = new ArrayList<Float>();
        suffix.add(lastCell.mass-lastCell.prevCells.get(i).mass);
        generateSequences(lastCell.prevCells.get(i), suffix, lastCell.score+lastCell.prevCells.get(i).score, scoreCutoff, results);
      }
    }
    return results;
  }
  
  
  private static void generateSequences(Cell current, ArrayList<Float> suffix, float cumScore, float scoreCutoff, ArrayList<ArrayList<Float>> results) {
    
    // base case
    if (current.scores.size()==0) {
      Collections.reverse(suffix);
      results.add(suffix);
      return;
    }
    
    for (int i=0; i<current.scores.size(); i++) {
      if (current.prevCells.get(i).bestScore+cumScore>=scoreCutoff) {
        //System.out.println("Recursive Score: " + current.prevCells.get(i).bestScore+cumScore);
        ArrayList<Float> extSuffix = new ArrayList<Float>(suffix);
        extSuffix.add(current.mass - current.prevCells.get(i).mass);
        generateSequences(current.prevCells.get(i), extSuffix, cumScore+current.prevCells.get(i).score, scoreCutoff, results);
      }
    }  
  }
  
}
