/**
 * 
 */
package msutil;

/**
 * The idea of this interface is that an implementing class can take an
 * Spectrum and spit out a (deep) copy of Spectrum with some properties
 * modified. Filters, recalibration and normalization are examples of
 * classes that should implement this interface. 
 * @author jung
 *
 */
public interface Reshape {
  
  /**
   * Apply this reshaping method for this Spectrum. 
   * @param s the spectrum to apply this operation
   * @return the new spectrum after applying the reshaping method. The input
   * spectrum is not changed.
   */
  public Spectrum apply(Spectrum s);

}
