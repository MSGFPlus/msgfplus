package suffixgraph.misc;

import java.util.Arrays;


/**
 * Space efficient implementation of a boolean/bit array.
 * @author jung
 *
 */
public class BitArray {
  int setItems;
  private int[] data;
  public BitArray() {
    data = new int[10];
    setItems =0;
  }
  
  private void increaseCapacity(int capacity) {
    int scaledCapacity = (capacity >>> 5) + 1;
    if (scaledCapacity >= data.length) {
      // expand array
      data = Arrays.copyOf(data, Math.max(scaledCapacity, data.length<<1));
    }
    return;  
  }
  
  /*
  public void unset(int index) {
    if (index >= size()) increaseCapacity(index);
    int scaledIndex = index >>> 5;   // index / 32 
    int offset = index & 0x1F;       // index mod 32}
    if (get(index))
      setItems--;
    data[scaledIndex] &= ~(1 << offset);
  }
  */
  
  public void set(int index) {
    if (index >= size()) increaseCapacity(index);
    int scaledIndex = index >>> 5;   // index / 32 
    int offset = index & 0x1F;       // index mod 32}
    if (!get(index))
      setItems++;
    data[scaledIndex] |= (1 << offset);
  }
  
  public boolean get(int index) {
    int scaledIndex = index >>> 5;   // index / 32 
    if (scaledIndex>=data.length) return false;
    int offset = index & 0x1F;       // index mod 32
    return ((data[scaledIndex] >>> offset) & 1) == 1;
  }
  
  public int size() {
    return data.length << 5;
  }
  
  public int getSetItems() {
    return setItems;
  }
}