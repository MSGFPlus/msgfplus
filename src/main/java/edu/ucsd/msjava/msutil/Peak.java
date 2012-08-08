package edu.ucsd.msjava.msutil;

import java.util.*;

/**
 * Representation of a peak in a spectrum object.
 * @author Sangtae Kim
 *
 */
public class Peak implements Comparable<Peak> {

  // required fields
  private int charge = 1;
  private float mz;
  private float intensity;
  
  // optional fields
  private int index = -1;
  private int rank = 151;
  
  
  /**
   * Constructor.
   * @param mz                         the m/z.
   * @param intensity                  the absolute intensity of this peak. 
   * @param charge                     the charge of this peak
   */
  public Peak(float mz, float intensity, int charge) {
    this.mz = mz;
    this.intensity = intensity;
    this.charge = charge;
//    if (charge!=0) this.charge = charge;	// commented out by Sangtae
  }
  
  
  /**
   * Gets the index of this peak.
   * @return the index of this peak of -1 if not initialized.
   */
  public int getIndex()                { return index; }

  
  /**
   * Gets the mz of this peak (as read from the file). 
   * @return the mass in Daltons.
   */
  public float getMz()               { return mz; }

  
  /**
   * Return the de-charged mass. 
   * @return (m/z - H) * charge.
   */
  public float getMass() {
    return  (mz - (float)Composition.H) * (float)charge; 
  }
  
  
  /**
   * Gets the intensity of this peak.
   * @return the intensity of this peak.
   */
  public float getIntensity()          { return intensity; }

  
  /**
   * Gets the charge of the peak.
   * @return the charge of the peak.
   */
  public int getCharge()               { return this.charge; }
  
  /**
   * Gets a new peak with different mass
   * @param modMass mass
   * @return a new peak with mass
   */
  public Peak getShiftedPeak(float mz)
  {
	  Peak newPeak = new Peak(mz, this.intensity, this.charge);
	  newPeak.rank = this.rank;
	  newPeak.index = this.index;
	  return newPeak;
  }
  
  /**
   * Sets the rank of this peak.
   * @param rank                       the rank of this peak.
   */
  public void setRank(int rank)        { this.rank = rank; }
  
  /**
   * Gets the rank of this peak.
   * @return the rank of this peaks or -1 if not initialized.
   */
  public int getRank()                 { return rank; }

  /**
   * Given the parent mass return the mass of the uncharged complement peak.
   * This assumes that the parent mass has no charge (H).
   * @param parentMass the deprotonated and decharged parent mass 
   * @return the deprotonated and decharged complement mass
   */
  public float getComplementMass(float parentMass) { 
    return parentMass - getMass(); 
  }

  
  /**
   * Sets the intensity of this peak.
   * @param intensity 
   */
  public void setIntensity(float intensity) {
    this.intensity = intensity;
  }

  
  /**
   * Sets the index of this peak.
   * @param index the index to set for this peak.
   */
  public void setIndex(int index) {
    this.index = index;
  }
  
  
  /**
   * Sets the mass of this peak to the given float.
   * @param mz the mass to set this peak to.
   */
  public void setMz(float mz) {
    this.mz = mz;
  }
  
  
  /**
   * Sets the charge of this peak.
   * @param charge the integer charge
   */
  public void setCharge(int charge) {
    this.charge = charge;
  }

  
  /**
   * Given ppm tolerance convert it to unit mass tolerance.
   * @param ppmTolerance the tolerance in ppm value
   * @return
   */
  public float toUnitTolerance(float ppmTolerance) {
    return getMass() * ppmTolerance / Constants.MILLION;
  }
  
  
  /**
   * 
   * @param shiftMass
   * @return
   */
  /*
  public Peak getShiftedPeak(float shiftMass) {
    Peak p = this.clone();
    p.mass += shiftMass;
    return p;
  }
  */

  
  
  /**
   * Compares this peak to another peak by mass. If the masses are equal, 
   * compare by intensity.
   */
  public int compareTo(Peak p) {
    if(mz > p.mz)                  return 1;
    if(p.mz > mz)                  return -1;
    
    if(intensity > p.intensity)        return 1;
    if(p.intensity > intensity)        return -1;
    
    return 0;
  }
  

  @Override
  public int hashCode()
  {
	 return (int)(mz+intensity+charge); 
  }
  
  /**
   * Checks the equality of this peak with another object.
   * @param obj the other object.
   * @return true if the intensities and masses are equal, false, otherwise.
   */
  @Override
  public boolean equals(Object obj)
  {
	  if(obj instanceof Peak)
		  return equals((Peak)obj);
	  return false;
  }
  
  /**
   * Checks the equality of this peak with another peak.
   * @param p the other peak.
   * @return true if the intensities and masses are equal, false, otherwise.
   */
  public boolean equals(Peak p) {
    // this might not be a good idea for floats
    return mz == p.mz && intensity == p.intensity && charge == p.charge;
  }
  
  
  /**
   * Calculates the absolute mass difference between 2 peaks. The m/z values 
   * are used for this method.
   * @param p1 the peak to subtract the mass from.
   * @param p2 the peak to subtract the mass by.
   * @return the mass difference in Daltons.
   */
  public static float getAbsoluteMassDiff(Peak p1, Peak p2) {
    return Math.abs(p1.mz - p2.mz);
  }

  
  /**
   * String representation of this peak. This is simply the mass followed by
   * its intensity.
   * @return mass, space, intensity string representation of this peak.
   */
  public String toString() {
    return mz + " " + intensity;
  }

  /**
   * Make a deep copy of this peak.
   * @return a deep copy of this peak
   */
  public Peak clone() {
    Peak p = new Peak(mz, intensity, charge);
    p.index = index;
    p.rank = rank;
    return p;
  }

  
  /**
   * Comparator to sort peaks by intensity. If the intensities are equal, sort
   * by mass
   * @author Sangtae Kim
   *
   */
  public static class IntensityComparator implements Comparator<Peak> {

    
    /**
     * Dictates ordering of peaks by intensity
     * @param p1 first peak to compare.
     * @param p2 second peak to compare.
     * @return 1 if p1 > p2, -1 if p2 > p1 and 0 if they are equal.
     */
    public int compare(Peak p1, Peak p2) {
      if(p1.intensity > p2.intensity)  return 1;
      if(p2.intensity > p1.intensity)  return -1; 
      
      if(p1.mz > p2.mz)            return 1;
      if(p2.mz > p1.mz)            return -1;
      
      return 0;
    }
    
    
    /**
     * Dictates equality of two peaks.
     * @param p1 first peak to compare.
     * @param p2 second peak to compare.
     * @return true, if both peaks have the same mass and intensity; false
     *         otherwise.
     */
    public boolean equals(Peak p1, Peak p2) {
      // this might not be a good idea because of float errors
      return p1.mz == p2.mz && p1.intensity == p2.intensity;
    }
  }
  
  /**
   * Dictates the order of peaks by mass.
   * @author Sangtae Kim
   *
   */
  public static class MassComparator implements Comparator<Peak> {

    /**
     * Comparison function by mass.
     * @param p1 first peak.
     * @param p2 second peak.
     * @return 1 if p1 > p2, -1 if p2 > p1 and 0 if they are equal.
     */
    public int compare(Peak p1, Peak p2) {
      return p1.compareTo(p2);
    }

    /**
     * Equality method.
     * @param p1 first peak.
     * @param p2 second peak.
     * @return true if their masses and intensities are equal, false otherwise.
     */
    public boolean equals(Peak p1, Peak p2) {
      return p1.equals(p2);
    }

  }
  
  
  
  /**
   * Creates a new peak with the same parameters as the current peak, but with
   * a mass offset given.
   * @param offset the offset to add
   * @return a peak object such that the getMass methods this return this.getMass()+offset
   *         as mass.
   */
  public Peak duplicate(float offset) {
    float mzOffset = offset / this.charge;
    return new Peak(mz+mzOffset, this.intensity, this.charge);
  }
  
}





