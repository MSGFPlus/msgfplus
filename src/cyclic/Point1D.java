package cyclic;

/**
 * This is the basic point data structure representing a spectral peak.
 * @author jung
 *
 */
public class Point1D implements Comparable<Point1D> {

  private static int count = 0;
  
  private float value;
  private float weight;
  private int index;
  private int id;
  
  /**
   * Constructor taking the required parameters for the Point1D object.
   * @param value the value (usually the mass of the point)
   * @param weight the weight for this item (usually the intensity or score of the peak)
   * @param index the index of the peak in the original spectrum . Use this to store any extra information.
   */
  public Point1D(float value, float weight, int index) {
    this.value = value;
    this.weight = weight;
    this.index = index;
    this.id = count++;
  }
  
  /**
   * Getter method.
   * @return the value of this point
   */
  public float getValue() { return this.value; }
  
  /**
   * Getter method.
   * @return the weight of this point
   */
  public float getWeight() { return this.weight; }
  
  /**
   * Getter method.
   * @return the peak index (or other information stored by this point
   */
  public int getIndex() { return this.index; }
  
  @Override
  public int compareTo(Point1D o) {
    if (this.value > o.value) return 1;
    if (this.value < o.value) return -1;
    
    if (this.id > o.id) return 1;
    if (this.id < o.id) return -1;
    
    return 0;
  }

}
