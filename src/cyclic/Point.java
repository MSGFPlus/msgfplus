package cyclic;

import java.util.Comparator;

import msutil.Spectrum;

/**
 * This represents an N-dimensional point in space.
 * @author jung
 *
 */
public class Point implements Comparable<Point> {

  /**
   * Helper method to determine what is the best rotation for this sequence of 
   * masses. The best rotation is always the one lexicographically minimal.
   * @param masses the ArrayList of masses
   * @return the index of the minimal rotation
   */
  private static int getMinRotation(float[] masses) {
    int minRot = 0;
    for (int i=1; i<masses.length; i++) {
      for (int subIndex=0; subIndex<masses.length; subIndex++) {
        float mass1 = masses[(subIndex+i)%masses.length];
        float mass2 = masses[(subIndex+minRot)%masses.length];
        if (mass1<mass2) {
          minRot = i;
          break;
        }
        else if (mass2<mass1) {
          // minRot = minRot;
          break;
        }
      }
    }
    return minRot;
  }
  
  
  private float[] masses;
  private int[] indices;
  private int rotation;
  private float weight;
  //private Spectrum spec;
  
  public Point(Spectrum spec, int[] indices, float weight) {
    this.masses = new float[indices.length];
    int i = 1;
    for (i=1; i<indices.length; i++) {
      this.masses[i-1] = spec.get(indices[i]).getMass() - spec.get(indices[i-1]).getMass();
    }
    this.masses[i-1] = spec.get(indices[0]).getMass() - spec.get(indices[i-1]).getMass() + spec.getParentMass();  
    this.indices = indices;
    this.rotation = getMinRotation(this.masses);
    this.weight = weight;
    //this.spec = spec;
  }
  
  public Point(Spectrum spec, int[] indices, float[] masses, float weight) {
    this.masses = masses;
    this.indices = indices;
    this.rotation = getMinRotation(this.masses);
    this.weight = weight;
    //this.spec = spec;
  }
  
  public int size() {
    return this.indices.length;
  }
  
  public float getMassAt(int index) {
    return masses[(index+this.rotation)%this.size()];
  }
  
  public float getWeight() { return this.weight; }
  
  
  /**
   * Calculate the euclidean distance of this point to the other point.
   * @param o the other Point object of the same dimension.
   * @return the distance between these 2 points
   */
  public float getDistance(Point o) {
    float cumDistance = 0.0f;
    for (int i=0; i<this.masses.length; i++) {
      float diff = this.getMassAt(i)-o.getMassAt(i);
      cumDistance = diff*diff;
    }
    return (float)Math.sqrt(cumDistance);
  }
  
  /**
   * Calculate the euclidean distance of this point to the other point. This is
   * an optimized method that will immediately return if the cumulative distance
   * is greater than maxDistance.
   * @param o the other Point object of the same dimension.
   * @param maxDistance the maximum distance  
   * @return the distance between these points if less than or equal to maxDistance, -1 otherwise.
   */
  public float getDistance(Point o, float maxDistance) {
    float cumDistance = 0.0f;
    float distSq = maxDistance*maxDistance;
    for (int i=0; i<this.masses.length; i++) {
      float diff = this.getMassAt(i)-o.getMassAt(i);
      cumDistance = diff*diff;
      if (cumDistance > distSq) return -1;
    }
    return (float)Math.sqrt(cumDistance);
  }

  /**
   * Return the maximum distance of these 2 points for all dimensions.
   * @param o the other Point object of the same dimension.
   * @return the distance
   */
  public float getMaxDistance(Point o) {
    float maxDistance = 0.0f;
    for (int i=0; i<this.masses.length; i++) {
      float diff = Math.abs(this.getMassAt(i)-o.getMassAt(i));
      if (maxDistance < diff) {
        maxDistance = diff;
      }
    }
    return maxDistance;
  }
  
  /**
   * Return the maximum distance of these 2 points for all dimensions.
   * @param o the other Point object of the same dimension.
   * @param maxDistance the maximum distance.
   * @return the distance or -1 if it exceeds maxDistance
   */
  public float getMaxDistance(Point o, float maxDistance) {
    float max = 0.0f;
    for (int i=0; i<this.masses.length; i++) {
      float diff = Math.abs(this.getMassAt(i)-o.getMassAt(i));
      if (max < diff) {
        max = diff;
      }
      if (max > maxDistance) return -1;
    }
    return maxDistance;
  }
  
  /**
   * Return the first peak index accounting for the rotation.
   * @return the first peak index.
   */
  public int getFirstIndex() {
    return this.indices[this.rotation%this.size()];  
  }
  
  @Override
  public int compareTo(Point o) {
    for (int i=0; i<size(); i++) {
      if (this.getMassAt(i)>o.getMassAt(i)) return 1;
      if (this.getMassAt(i)<o.getMassAt(i)) return -1;
    }
    return 0;
  }
  
  @Override
  public String toString() {
    StringBuffer sb = new StringBuffer();
    for (int i=0; i<this.size(); i++) {
      sb.append(String.format("%.2f ", this.getMassAt(i)));
    }
    return sb.toString();
  }
  
  /**
   * Sorting by weight in descending order 
   * @author jung
   *
   */
  private static class PointWeightComparator implements Comparator<Point> {
    @Override
    public int compare(Point a, Point b) {
      if (a.weight > b.weight) return -1;
      if (a.weight < b.weight) return 1;
      return 0;
    }
  }
  public static final PointWeightComparator pointWeightComparator = new PointWeightComparator();
  
}
