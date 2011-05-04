package suffixtree.edges;

import java.util.Arrays;


/**
 * This Object is to create edges for querying a graph.
 * @author jung
 *
 */
public class ByteEdge extends Edge {

  int label;
  
  /**
   * Constructor taking a byte array.
   * @param sequence the byte array as the label of this edge.
   */
  public ByteEdge(int label) {
    this.label = label;
  }
  
  
  public ByteEdge(byte[] label) {
    if (label.length>4) 
      System.err.println("Attempting to create a QueryEdge with a label greater than 4 bytes... truncating.");
    this.label = 0;
    Arrays.sort(label);
    for (int i=0; i<label.length; i++) {
      this.label = this.label<<8;
      this.label |= label[i];
    }
  }
  
  @Override
  public int getLabelAt(int offset) {
    return label;
  }

  @Override
  public int size() {
    return 1;  
  }
  
  @Override
  public int length() {
    int mask = 0xFF;
    for (int i=0; i<4; i++) {
      if ((label&mask)==0) return i;
      mask = mask<<8;
    }
    return 4;
  }

  @Override
  public Edge split(int offset) {
    System.err.println("Unsupported operation: Split query edge");
    return null;
  }

  @Override
  public int getEnd() {
    System.err.println("Unsupported operation");
    System.exit(-1);
    return 0;
  }

  @Override
  public int getStart() {
    System.err.println("Unsupported operation");
    System.exit(-1);
    return 0;
  }
  
  
  /**
   * Expand the compressed label of this item to a byte array
   * @return the byte array with individual bytes decompressed.
   */
  public byte[] toByteArray() {
    if ((0xFFFFFF00&this.label)==0) {
      byte[] result =  {(byte)this.label};
      return result;
    }
    
    if ((0xFFFF0000&this.label)==0) {
      byte[] result =  {(byte)(this.label>>>8), (byte)this.label};
      return result;
    }
    
    if ((0xFF000000&this.label)==0) {
      byte[] result =  {(byte)(this.label>>>16), (byte)(this.label>>>8), (byte)this.label};
      return result;
    }
    
    byte[] result =  {(byte)(this.label>>>24), (byte)(this.label>>>16), (byte)(this.label>>>8), (byte)this.label};
    return result;
  }
  
  
  /**
   * Remove the given bytes from the current edge and return the edge after
   * removing the given bytes.
   * @param bytes the bytes to remove. The length of this array must be strictly 
   *              less than the length of this edge
   * @return a composite edge resulting from removing the bytes from the current
   *         edge. null is returned in case that bytes is not contained in the
   *         current edge.
   */
  public ByteEdge remove(byte[] bytes) {
    
    Arrays.sort(bytes);
    int foundCount = 0;
    
    int mask = 0xFF000000;
    int result = 0;
    // iterate over the elements of the current edge
    for (int i=0; i<32; i+=8) {
      byte item = (byte)((this.label&mask)>>>(24-i));
      if (item != 0) {
        
        if (foundCount<bytes.length && bytes[foundCount]==item) {
          foundCount++;          
        }
        else {
          // keep this item
          result = result<<8;
          result = result | item;
        }
      }
      mask = mask>>>8;
    }
    
    if (foundCount<bytes.length) {
      // some items do not exist in the current edge
      return null;
    }
    
    return new ByteEdge(result);
  }


  @Override
  public int mass() {
    return label;
  }
}
