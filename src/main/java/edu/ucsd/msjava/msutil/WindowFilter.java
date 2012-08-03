package edu.ucsd.msjava.msutil;


/**
 * This filtering method guarantees that for any window of the determined size
 * placed around a given peak, the peak is ranked better than a determined 
 * parameter.
 * @author jung
 *
 */
public class WindowFilter implements Reshape {
  
  // the number of peaks to retain.
  private int top;
  
  // the width of the window. +/- Daltons.
  private float window;
  
  
  
  /**
   * Constructor.
   * @param top the number of peaks to keep per window. 1-based rank.
   * @param window the window size in Daltons. The window will be taken to the
   *               left and right of the amount in Daltons of the specified parameter
   */
  public WindowFilter(int top, float window) {
    this.top = top;
    this.window = window;
  }

  /**
   * Getter method.
   * @return the number of top peaks for this window.
   */
  public int getTop() { return top; }

  /**
   * Getter methods.
   * @return the window width used to create this window filter.
   */
  public float getWindow() { return window; }
  
  public Spectrum apply(Spectrum s) {
    
    // select each peak if it is top n within window (-window,+window) around it
    Spectrum retSpec = (Spectrum)s.clone();
    retSpec.clear();    // remove all peaks 
   
    for(int peakIndex = 0; peakIndex < s.size(); peakIndex++) {
      int rank = 1;
      
      Peak thisPeak = s.get(peakIndex);
      float thisMass = thisPeak.getMass();
      float thisInten = thisPeak.getIntensity();
      
      // move left
      int prevIndex = peakIndex-1;
      while(prevIndex >= 0) {
        Peak prevPeak = s.get(prevIndex);
        if(thisMass - prevPeak.getMass() > this.window)    break;
        if(prevPeak.getIntensity() > thisInten)            rank++;
        prevIndex--;
      }

      // move right
      int nextIndex = peakIndex+1;
      while(nextIndex < s.size()) {
        Peak nextPeak = s.get(nextIndex);
        if(nextPeak.getMass() - thisMass > this.window)    break;
        if(nextPeak.getIntensity() > thisInten)            rank++;
        nextIndex++;
      }
    
      if(rank <= this.top)   retSpec.add(thisPeak);
    }
    return retSpec;
  }
}
