package edu.ucsd.msjava.msutil;


/**
 * Retain a fixed number of peaks ranked by intensity.
 * @author jung
 *
 */
public class RankFilter implements Reshape {
  
  // the number of peaks to retain.
  private int top;
  

  
  /**
   * Constructor.
   * @param top the number of peaks to keep in the 1-based rank.
   */
  public RankFilter(int top) {
    this.top = top;
  }
  
  
  
  /**
   * Reshape the given spectrum by discarding all peaks that are below a given
   * rank.
   */
  public Spectrum apply(Spectrum s) {
    
    // select each peak if it is top n within window (-window,+window) around it
    Spectrum retSpec = (Spectrum)s.clone();
    s.setRanksOfPeaks();
    retSpec.clear();    // remove all peaks 
   
    for(Peak thisPeak : s) {
      if(thisPeak.getRank() <= this.top)   retSpec.add(thisPeak.clone());
    }
    
    return retSpec;
  }

}
