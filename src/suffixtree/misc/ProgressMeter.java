package suffixtree.misc;

import java.io.PrintStream;


/**
 * This class prints out the progress meter for the building big datastructures
 * @author jung
 *
 */
public class ProgressMeter {
  
  private final static int STEP = 5;
  private PrintStream out;
  private long maxCount;
  private float nextCount;
  private float step;
  
  public ProgressMeter(String message, long maxCount, PrintStream out) {
    this.out = out;
    this.out.printf("%s 0%%", message);
    this.maxCount = maxCount;
    this.step = maxCount * STEP / 100.0f; 
    this.nextCount = this.step;
  }
  
  public void update(long count) {
    if (count > this.nextCount) {
      this.out.printf(" %d%%", (int)Math.ceil(count*100/this.maxCount));
      this.nextCount += this.step;
    }
  }
  
  public void update(long count, String message) {
    if (count > this.nextCount) {
      this.out.printf(" %d%% %s", (int)Math.ceil(count*100/this.maxCount), message);
      this.nextCount += this.step;
    }
  }
  
}
