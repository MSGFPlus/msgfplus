package cyclic;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

import msutil.Spectrum;

public class Cluster1D implements Comparable<Cluster1D> {
  
  /**
   * Holds 2 clusters, A and B.
   * @author jung
   *
   */
  private static class Pair implements Comparable<Pair> {

    private static int pairCount = 0;
    private Cluster1D a;
    private Cluster1D b;
    private float distance;
    private int id;
    
    public Pair(Cluster1D a, Cluster1D b, float distance) {
      this.a = a;
      this.b = b;
      this.distance = distance;
      this.id = pairCount++;
    }
    
    @Override
    public int compareTo(Pair o) {
      if (this.distance > o.distance) return 1;
      if (this.distance < o.distance) return -1;
      return this.id - o.id;
    }
  }
  
  private static int count = 0;
  
  private ArrayList<Point1D> points;
  private float center;
  private float weight;
  private int id;

  public Cluster1D(Point1D p) {
    this.center = p.getValue();
    this.weight = p.getWeight();
    this.points = new ArrayList<Point1D>();
    this.points.add(p);
    this.id = count++;
  }
  
  private Cluster1D(float center, float weight, ArrayList<Point1D> points) {
    this.center = center;
    this.weight = weight;
    this.points = points;
    this.id = count++;
  }
  
  public int size() {
    return this.points.size();
  }
  
  public ArrayList<Point1D> getPoints() { return this.points; }
  
  public float getCenter() { return this.center; }
  
  public float getWeight() { return this.weight; }
  
  public Cluster1D merge(Cluster1D other) {
    float newWeight = this.weight + other.weight;
    this.center = (this.center*this.weight + other.center*other.weight) / newWeight;
    this.weight = newWeight;
    this.points.addAll(other.points);
    return new Cluster1D(this.center, this.weight, this.points);
  }
  
  @Override
  public String toString() {
    return String.format("Cluster1D [%.2f %.2f]", this.center, this.weight);  
  }
  
  
  @Override
  public int compareTo(Cluster1D o) {
    if (this.center > o.center) return 1;
    if (this.center < o.center) return -1;
    
    if (this.id > o.id) return 1;
    if (this.id < o.id) return -1;
    
    return 0;
  }
  
  public static Set<Cluster1D> make1DClusters(Spectrum spec, float minMass, float tolerance, int minClusterSize) {
    return cluster(convolution(spec, minMass), tolerance, minClusterSize);
  }
  
  private static ArrayList<Point1D> convolution(Spectrum spec, float minMass) {
    ArrayList<Point1D> points = new ArrayList<Point1D>();
    for (int i=0; i<spec.size(); i++) {
      for (int j=i+1; j<spec.size(); j++) {
        float diff1 = spec.get(j).getMass() - spec.get(i).getMass();
        float diff2 = spec.getParentMass()-diff1;
        float score = spec.get(j).getIntensity() + spec.get(i).getIntensity();
        if (diff1 < diff2) {
          points.add(new Point1D(diff1, score, 0));      
        }
        else {
          points.add(new Point1D(diff2, score, 0));
        }
      }
    }
    return points;
  }
  
  /**
   * This is the main (greedy) clustering routine for 1D clustering.
   * @param points the list of points to clusters
   * @param maxDistance the maximum radius for the clustering 
   * @param minClusterSize the minimum size of the cluster
   * @return the Set of resulting clusters
   */
  public static Set<Cluster1D> cluster(ArrayList<Point1D> points, float maxDistance, int minClusterSize) {
    
    TreeSet<Pair> sortedDists = new TreeSet<Pair>();
    
    TreeSet<Cluster1D> clusters = new TreeSet<Cluster1D>();
    for (Point1D p : points) {
      //System.out.println(p.getValue());
      clusters.add(new Cluster1D(p));
    }
    
    Iterator<Cluster1D> it = clusters.iterator();
    
    Cluster1D prevCluster = it.next();
    while (it.hasNext()) {
      Cluster1D currentCluster = it.next();
      if (currentCluster.center-prevCluster.center<=maxDistance) {
        sortedDists.add(new Pair(prevCluster, currentCluster, currentCluster.center-prevCluster.center));
      }
      prevCluster = currentCluster;
    }
    
    while (sortedDists.size() > 0) {
      Pair target = sortedDists.pollFirst();
      while (!clusters.contains(target.a) || !clusters.contains(target.b)) {
        if (sortedDists.size() > 0) {
          // this is a valid pair
          target = sortedDists.pollFirst();
        }
        else {
          it = clusters.iterator();
          while (it.hasNext()) if (it.next().size()<minClusterSize) it.remove();
          return clusters;
        }
      }
      
      // remove both clusters from the final list
      clusters.remove(target.a);
      clusters.remove(target.b);

      // create and add the merged cluster to the distances
      Cluster1D merged = target.a.merge(target.b);
      Cluster1D smaller = clusters.floor(merged);
      Pair minPair = null;
      if (smaller!=null) {
        minPair = new Pair(smaller, merged, merged.center-smaller.center);
      }
      Cluster1D larger = clusters.ceiling(merged);
      if (larger!=null) {
        Pair temp = new Pair(merged, larger, larger.center-merged.center);
        if (minPair==null || temp.compareTo(minPair)<0) minPair = temp;
      }
      if (minPair==null) {
        // we had no clusters
        clusters.add(merged);
        return clusters;
      }
      
      clusters.add(merged);
      if (minPair.distance<=maxDistance) {
        sortedDists.add(minPair);
      }
    }
    
    it = clusters.iterator();
    while (it.hasNext()) if (it.next().size()<minClusterSize) it.remove();
    return clusters;
  }

  
  public static Cluster1DWeightComparator weightComparator = new Cluster1DWeightComparator();
  /**
   * Allow sort by weight of the clusters. Higher weights come first
   * @author jung
   *
   */
  private static class Cluster1DWeightComparator implements Comparator<Cluster1D> {

    @Override
    public int compare(Cluster1D a, Cluster1D b) {
      if (a.weight > b.weight) return -1;
      if (a.weight < b.weight) return 1;
      
      if (a.id > b.id) return -1;
      if (a.id < b.id) return 1;
      return 0;
    }
    
  }
}
