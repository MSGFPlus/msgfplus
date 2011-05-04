package suffixgraph.misc;


/**
 * This is a space efficient implementation of a priority queue. The backing
 * data structure is an array of primitive types long. This is implemented as
 * a Min Heap. This allows for duplicate elements.
 * @author jung
 *
 */
public class LongPriorityQueue {

  private int size;
  private int levels;   // the number of levels in this heap
  private long[] items;
  
  public LongPriorityQueue() {
    this.levels = 8;
    this.size = 0;
    this.items = new long[(1<<levels)-1];
    for (int i=0; i<items.length; items[i++] = Long.MAX_VALUE);
  }
  
  
  /**
   * Remove minimum from the heap.
   * @return the minimum element. If the queue is empty it will return Long.MAX_VALUE.
   */
  public long poll() {
    long value = items[0];
    items[0] = items[size-1];
    items[size-1] = Long.MAX_VALUE;
    size--;
    tickleDown(0);
    return value;
  }
  
  
  /**
   * Add a given element into the queue.
   * @param e the element to add.
   */
  public void add(long e) {
    // check that the array can fit this item, the duplication will guarantee
    // that the heap is a complete tree with maximum elements filling the 
    // unused leafs
    if (size >= items.length) {
      long[] temp = new long[(1<<(levels+1))-1];   // double the array size
      System.arraycopy(items, 0, temp, 0, items.length);
      items = temp;
      for (int i=1<<levels; i<items.length; items[i++] = Long.MAX_VALUE);
      levels++;
    }
    
    items[size] = e;
    bubbleUp(size);
    this.size++;
  }

  
  /**
   * Return the number of elements in this data structure.
   * @return the number of elements in this object.
   */
  public int size() {
    return this.size;
  }
  
  
  /**
   * Check whether there are elements in this object.
   * @return false if there size > 0, true otherwise.
   */
  public boolean isEmpty() {
    return size==0;
  }
  
  
/***** HELPING METHODS *****/
  /**
   * Simply loop to maintain the properties of the minimum heap.
   * @param targetIndex the index of the newly inserted item.
   */
  private void bubbleUp(int targetIndex) {
    while (targetIndex>0) {
      int parentIndex = (targetIndex-1)>>>1;
      if (items[parentIndex] > items[targetIndex]) {
        // do a swap
        long temp = items[targetIndex];
        items[targetIndex] = items[parentIndex];
        items[parentIndex] = temp;
      }
      else {
        // we are done
        break;
      }
      targetIndex = parentIndex;
    }
  }
  
  /**
   * Helper method that maintains the min heap properties after a remove min.
   * @param targetIndex the index of the remove min.
   */
  private void tickleDown(int targetIndex) {
    
    int rChildIndex = (targetIndex+1)<<1;
    int lChildIndex = rChildIndex-1;
    while (rChildIndex<items.length) {

      int minIndex = 0;
      if (items[lChildIndex]<items[rChildIndex]) { 
        minIndex = lChildIndex;
      }
      else {
        minIndex = rChildIndex;
      }
      
      if (items[targetIndex]>items[minIndex]) {
        // do a swap
        long temp = items[targetIndex];
        items[targetIndex] = items[minIndex];
        items[minIndex] = temp;
        targetIndex = minIndex;
      }
      else {
        break;
      }
      
      // update the pointers
      rChildIndex = (targetIndex+1)<<1;
      lChildIndex = rChildIndex-1;
    }
  }
  
  
}
