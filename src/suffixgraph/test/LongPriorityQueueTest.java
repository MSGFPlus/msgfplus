package suffixgraph.test;

import java.util.Arrays;
import java.util.PriorityQueue;
import java.util.Random;

import suffixgraph.misc.LongPriorityQueue;

public class LongPriorityQueueTest {

  /**
   * The testing method against a regular PriorityQueue.
   * @param args
   */
  public static void main(String[] args) {
    int maxElements = 10000000;   // number of elements to test
    
    Random gen = new Random();
    
    long[] elements = new long[maxElements];
    for (int i=0; i<maxElements; elements[i++]=gen.nextLong());
    
    PriorityQueue<Long> otherQ = new PriorityQueue<Long>();
    long startTime = System.currentTimeMillis();
    for (long nextLong : elements) {
      otherQ.add(nextLong);
    }
    System.out.printf("Time to insert %d elements on standard queue %.2fs\n", maxElements, (System.currentTimeMillis()-startTime)/1000.0);
    
    LongPriorityQueue myQ = new LongPriorityQueue();
    startTime = System.currentTimeMillis();
    for (long nextLong : elements) {
      myQ.add(nextLong);
    }
    System.out.printf("Time to insert %d elements on custom queue %.2fs\n", maxElements, (System.currentTimeMillis()-startTime)/1000.0);
       
    Arrays.sort(elements);
    startTime = System.currentTimeMillis();
    for (int i=0; i<maxElements; i++) {
      if (otherQ.poll()!=elements[i]) { 
        System.out.printf("Standard min queue error");
        System.exit(-1);
      }
    }
    System.out.printf("Time to remove %d elements on standard queue %.2fs\n", maxElements, (System.currentTimeMillis()-startTime)/1000.0);
    
    startTime = System.currentTimeMillis();
    for (int i=0; i<maxElements; i++) {
      if (myQ.poll()!=elements[i]) { 
        System.out.printf("Custom min queue error");
        System.exit(-1);
      }
    }
    System.out.printf("Time to remove %d elements on custom queue %.2fs\n", maxElements, (System.currentTimeMillis()-startTime)/1000.0);

  }
  
}
