package msutil;


/**
 * This is root class for anything that has a mass. 
 * @author jung
 *
 */
public abstract class Matter implements Comparable<Matter> {
	
   /**
   * Get the mass of this object.
   * @return the mass in Daltons of this object.
   */
  public abstract float getMass();
  
  /**
   * Get the accurate (double-precision) mass of this object
   * @return
   */
  public double getAccurateMass()
  {
	  return getMass();
  }
  
  /**
   * Get the nominal (integer) mass of this object.
   * @return the nominal mass in Daltons of this object
   */
  public abstract int getNominalMass();
  
  /**
   * Defines the ordering of amino acids. Order by mass.
   * @param other massive to compared it to.
   * @return 1 if this is greater than the other, -1 if the other is greater 
   *         than this and 0 if they are equal.
   */
  public int compareTo(Matter other) {
    if(this.getMass() > other.getMass())         return 1;
    if(other.getMass() > this.getMass())         return -1;
    return 0;
  }
	  
  /**
   * Standard string representation of this object.
   * @return mass in string with 2 significant figures.
   */
  public String toString() {
    return String.format("[%.2f]", getMass()); 
  }
  
  public abstract boolean equals(Object obj);
}
