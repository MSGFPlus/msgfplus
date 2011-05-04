package cyclic;

public class Constants {
  
  // the maximum number of points to take for a cluster
  public static final int MAX_CLUSTER_PEAKS = 20; 
  
  // the minimum number of members in a cluster
  public static final int MIN_CLUSTER_SIZE = 4; 
  
  // the minimum peak count before introducing a consensus peak
  public static final int MIN_PEAK_COUNT = 4;
  
  // the minimum separation for adjacent valid peaks
  public static final float MIN_DISTANCE = 57;
  
  // the upper limit for single masses
  public static final float MAX_DISTANCE = 200;
  
  // the number of top single convolution masses to use
  public static final int MAX_CONV_MASSES = 50;
}
