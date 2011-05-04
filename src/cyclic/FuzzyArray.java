package cyclic;

import java.util.Arrays;

public class FuzzyArray {

  /**
   * Emulates a binary search but within a certain tolerance.
   * @param array the float array to do binary search on
   * @param key the key to look for
   * @param tolerance the maximum tolerance to the match
   * @return the index of the closest item within tolerance or a negative number if no match.
   */
  public static int search(float[] array, float key, float tolerance) {
    
    int matchIndex = Arrays.binarySearch(array, key);
    if (matchIndex>=0) return matchIndex;
    
    // convert the match index into a positive number
    matchIndex = -matchIndex-1;

    // check the distance to the item on the right
    float rightDelta = Float.MAX_VALUE;
    if (matchIndex<array.length) {
      rightDelta = Math.abs(array[matchIndex]-key);
    }
    
    // check the distance to the item on the left
    float leftDelta = Float.MAX_VALUE;
    if (matchIndex>0) {
      leftDelta = Math.abs(key-array[matchIndex-1]);
    }  
    
    if (rightDelta<leftDelta) {
      if (rightDelta<=tolerance) return matchIndex;
      else                       return -matchIndex-1;
    }
    
    else {
      if (leftDelta<=tolerance) return matchIndex-1; 
      else                      return -(matchIndex-1)-1;
    }
  }
  
  
  /**
   * Debugging routine.
   * @param args
   */
  public static void main(String[] args) {
    float[] array = {0, 1, 2, 3, 4, 5, 6};
    System.out.println(search(array, 5.05f, 0.1f));
  }
  
}
