package cyclic;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;
import java.util.TreeSet;

import msutil.Spectrum;

/**
 * Holds a series of points
 * @author jung
 *
 */
public class Cluster implements Comparable<Cluster> {
  
  /**
   * Holds 2 clusters, A and B.
   * @author jung
   *
   */
  private static class Pair implements Comparable<Pair> {

    private static int pairCount = 0;
    private Cluster a;
    private Cluster b;
    private float distance;
    private int id;
    
    /*
    public Pair(Cluster a, Cluster b) {
      this.a = a;
      this.b = b;
      this.distance = a.getDistance(b);
      this.id = pairCount++;
    }
    */
    
    public Pair(Cluster a, Cluster b, float distance) {
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
  
  private float[] centers;
  private float weight;
  private ArrayList<Point> points;
  private int id;
  
  private Cluster(float[] centers, float weight, ArrayList<Point> points) {
    this.centers = centers;
    this.weight = weight;
    this.points = points;
    this.id = count++;
  }
  
  public Cluster(Point p) {
    this.centers = new float[p.size()];
    for (int i=0; i<this.centers.length; i++) {
      this.centers[i] = p.getMassAt(i);
    }
    this.weight = p.getWeight();
    this.points = new ArrayList<Point>();
    this.points.add(p);
    this.id = count++;
  }
  
  public void add(Point p) {
    float newWeight = this.weight + p.getWeight();
    for (int i=0; i<centers.length; i++) {
      centers[i] = (centers[i]*weight + p.getMassAt(i)*p.getWeight()) / newWeight;
    }
    this.weight = newWeight;
    this.points.add(p);
    this.id = count++;
  }
  
  public int getId() {
    return this.id;  
  }
  
  /**
   * Getter method that returns the centers of this object
   * @return the centers 
   */
  public float[] getCenters() {
    return this.centers;  
  }
  
  /**
   * Get the point members from this object
   * @return the list of points.
   */
  public ArrayList<Point> getPoints() {
    return this.points;
  }
  
  /**
   * Calculate the euclidean distance of this point to the other point.
   * @param o the other Point object of the same dimension.
   * @return
   */
  public float getDistance(Cluster o) {
    float cumDistance = 0.0f;
    for (int i=0; i<this.centers.length; i++) {
      float diff = this.centers[i]-o.centers[i];
      cumDistance = diff*diff;
    }
    return (float)Math.sqrt(cumDistance);
  }

  /**
   * Calculate the euclidean distance of this cluster to the other cluster. This is
   * an optimized method that will immediately return if the cumulative distance
   * is greater than maxDistance.
   * @param o the other Cluster object of the same dimension.
   * @param maxDistance the maximum distance  
   * @return the distance between these clusters if less than or equal to maxDistance, -1 otherwise.
   */
  public float getDistance(Cluster o, float maxDistance) {
    float cumDistance = 0.0f;
    float distSq = maxDistance*maxDistance;
    for (int i=0; i<this.centers.length; i++) {
      float diff = this.centers[i]-o.centers[i];
      cumDistance = diff*diff;
      if (cumDistance > distSq) return -1;
    }
    return (float)Math.sqrt(cumDistance);
  }
  
  /**
   * Return the maximum distance of these 2 points for all dimensions.
   * @param o the other Point object of the same dimension.
   * @return
   */
  public float getMaxDistance(Cluster o) {
    float maxDistance = 0.0f;
    for (int i=0; i<this.centers.length; i++) {
      float diff = Math.abs(this.centers[i]-o.centers[i]);
      if (maxDistance < diff) {
        maxDistance = diff;
      }
    }
    return maxDistance;
  }
  
  /**
   * Add this cluster with the other cluster. The calling object is modified.
   * @param o the other cluster object
   */
  public void add(Cluster c) {
    float newWeight = this.weight + c.weight;
    for (int i=0; i<centers.length; i++) {
      centers[i] = (this.centers[i]*weight + c.centers[i]*c.weight) / newWeight;
    }
    this.weight = newWeight;
    this.points.addAll(c.points);
  }
  
  /**
   * Merge this cluster with the other cluster. The calling object is modified.
   * @param o the other cluster object
   */
  public Cluster merge(Cluster c) {
    float newWeight = this.weight + c.weight;
    for (int i=0; i<centers.length; i++) {
      centers[i] = (this.centers[i]*weight + c.centers[i]*c.weight) / newWeight;
    }
    this.weight = newWeight;
    this.points.addAll(c.points);
    return new Cluster(this.centers, this.weight, this.points);
  }
  
  @Override
  public int compareTo(Cluster o) {
    for (int i=0; i<centers.length; i++) {
      if (this.centers[i] > o.centers[i]) {
        return 1;
      }
      if (this.centers[i] < o.centers[i]) {
        return -1;
      }
    }
    return 0;
  }
  
  /**
   * Return the number of points that form this cluster.
   * @return the number of points in this cluster
   */
  public int size() {
    return this.points.size();
  }
  
  /**
   * Return the 2D space hash of this point. This is the grid id in which this
   * cluster would fall, given the step, which the width of the grid.
   * @param step the width of the grid to place the this cluster.
   * @return the long hash number width the upper 32 bits being the hash of the
   *         first coordinate and the lower 32 bits being the hash of the second 
   *         coordinate
   */
  public long get2DHash(float step) {
    return (((long)(this.centers[0]/step))<<32) | ((long)(this.centers[1]/step));
  }

  /**
   * Get the near hashes, which is rectangle of 9 cells.
   * @param step the width of the cell
   * @return an array of near hashes
   */
  public long[] getNear2DHashes(float step) {
    long x = (((long)(this.centers[0]/step))<<32);
    long y = ((long)(this.centers[1]/step));
    long ul = (x-1)|(y-1), uc = x|(y-1), ur = (x+1)|(y-1),
         ml = (x-1)|y, mc = x|y, mr = (x+1)|y,
         ll = (x-1)|(y+1), lc = x|(y+1), lr = (x+1)|(y+1);
    long[] ret = {ul, uc, ur, ml, mc, mr, ll, lc, lr};
    return ret;
  }
  
  /**
   * Reverse order comparator. Put the higher weight first.
   * @author jung
   *
   */
  public static class ClusterWeightComparator implements Comparator<Cluster> {
    @Override
    public int compare(Cluster a, Cluster b) {
      if (a.weight>b.weight) return -1;
      if (a.weight<b.weight) return 1;
      return 0;
    }
  }
  
  @Override
  public String toString() {
    StringBuffer sb = new StringBuffer();
    for (float center : this.centers) {
      sb.append(String.format("%.2f ", center));
    }
    sb.append("; ");
    sb.append(this.size());
    sb.append(" ; ");
    sb.append(this.weight);
    return sb.toString();
  }
  

  public static Set<Cluster> make2DClusters(Spectrum spec, float minMass, float tolerance, int minClusterSize) {
    return cluster(convolution2D(spec, minMass), tolerance, minClusterSize);
  }
  
  private static ArrayList<Point> convolution2D(Spectrum spec, float minMass) {
    
    // get the differences for the convolution
    ArrayList<Point> points = new ArrayList<Point>();
    for(int i=0; i<spec.size()-2; i++) {
      for(int j=i+1; j<spec.size()-1; j++) {
        
        float diff1 = spec.get(j).getMass() - spec.get(i).getMass();
        float score = spec.get(j).getIntensity() + spec.get(i).getIntensity();
        //System.out.println(spec.get(j).getMass() + " " + spec.get(j).getCharge());
        if (diff1 < minMass) continue; // don't care about small masses
        
        for (int k=j+1; k<spec.size(); k++) {
          
          float diff2 = spec.get(k).getMass() - spec.get(j).getMass();
          if (diff2 < minMass) continue; // mass is too small
          
          float diff3 = spec.get(i).getMass() - spec.get(k).getMass() + spec.getParentMass();
          if (diff3 < minMass) continue; // mass is too small

          int[] indices = {i, j ,k};
          float[] masses = {diff1, diff2, diff3};
          points.add(new Point(spec, indices, masses, score+spec.get(k).getIntensity()));
        }
      }
    }
    return points;
  }
  
  
  public static Set<Cluster> cluster(ArrayList<Point> points, float maxDistance, int minClusterSize) {
    TreeSet<Pair> sortedDists = new TreeSet<Pair>();
    
    // spread the points into a 2D hash
    HashMap<Long,Collection<Cluster>> spaceHash = new HashMap<Long,Collection<Cluster>>();
    HashSet<Cluster> clusters = new HashSet<Cluster>();
    
    for (Point p : points) {
      Cluster c = new Cluster(p);
      long key = c.get2DHash(maxDistance);
      if (!spaceHash.containsKey(key)) {
        spaceHash.put(key, new LinkedList<Cluster>());  
      }
      spaceHash.get(key).add(c);
      clusters.add(c);
    }
    //System.out.println("Done creating the spaceHash");  

    for (Cluster a : clusters) {
      for (long key : a.getNear2DHashes(maxDistance)) {
        if (spaceHash.containsKey(key)) {
          for (Cluster b : spaceHash.get(key)) {
            if (a.getId() < b.getId()) {
              float distance = a.getDistance(b, maxDistance);
              if (distance>=0) {
                Pair p = new Pair(a, b, distance);
                sortedDists.add(p);
              }
            }
          }
        }
      }
    }
    
    //System.out.println("Initial distances calculated... starting to cluster.");
    
    while (sortedDists.size() > 0) {
      Pair target = sortedDists.pollFirst();
      while (!clusters.contains(target.a) || !clusters.contains(target.b)) {
        if (sortedDists.size() > 0) {
          // this is a valid pair
          target = sortedDists.pollFirst();
        }
        else {
          Iterator<Cluster> it = clusters.iterator();
          while (it.hasNext()) if (it.next().size()<minClusterSize) it.remove();
          return clusters;
        }
      }
      
      clusters.remove(target.a);
      clusters.remove(target.b);
      Cluster merged = target.a.merge(target.b);
      
      // add the new merged cluster to the space hash
      long margedHash = merged.get2DHash(maxDistance);
      if (!spaceHash.containsKey(margedHash)) {
        spaceHash.put(margedHash, new LinkedList<Cluster>());  
      }
      spaceHash.get(margedHash).add(merged);
      clusters.add(merged);
      
      // add the new distances
      for (long key : merged.getNear2DHashes(maxDistance)) {
        if (spaceHash.containsKey(key)) {
          Iterator<Cluster> it = spaceHash.get(key).iterator();
          while (it.hasNext()) {
            Cluster other = it.next();
            if (!clusters.contains(other)) {
              // clean up the items in this cell because this item has been removed
              it.remove();
            }
            else {
              if (other.getId()!=merged.getId()) {
                float distance = merged.getDistance(other, maxDistance);
                if (distance>=0) {
                  sortedDists.add(new Pair(other, merged, distance));
                }
              }
            }
          }
        }
      }
      
      //System.out.println("PQ size: " + sortedDists.size() + ". Clusters size: " + clusters.size());
    }
    
    Iterator<Cluster> it = clusters.iterator();
    while (it.hasNext()) if (it.next().size()<minClusterSize) it.remove();
    return clusters;
  }
}
