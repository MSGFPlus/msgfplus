package suffixgraph.test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import msutil.AminoAcid;
import msutil.Composition;
import suffixgraph.nodes.Node;

public class NodeTest {

  /**
   * Print out statistics about memory usage and time to run ArrayList VS array.
   * @param iterations
   */
  private static void edgeInsertTimeMemoryTest(int iterations) {
    Random g = new Random();
    AminoAcid[] stdAA = AminoAcid.getStandardAminoAcids();
    
    long time = System.currentTimeMillis();
    ArrayList<Composition> keyOrder = new ArrayList<Composition>();
    ArrayList<Integer> nodes = new ArrayList<Integer>();
    for (int i=0; i<iterations; i++) {
      int comp = stdAA[(Math.abs(g.nextInt())%stdAA.length)].getComposition().getNumber();
      
      Composition keyComp = new Composition(comp);
      int insertIndex = Collections.binarySearch(keyOrder, keyComp);
      keyOrder.add((insertIndex<0) ? -insertIndex - 1 : insertIndex, keyComp);
      nodes.add(-3);
      
      if (i%100000 == 0) {
        long usedMem = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024);
        System.out.println("Used " + usedMem + " MB");
      }
    }
    System.out.println("-- Testing time for Collection method: " + (System.currentTimeMillis() - time)/1000.0 + "s");
    keyOrder = null;
    nodes = null;
    System.gc();
    
    time = System.currentTimeMillis();
    Node node = new Node();
    for (int i=0; i<iterations; i++) {
      int comp = stdAA[(Math.abs(g.nextInt())%stdAA.length)].getComposition().getNumber();
      node.addEdge(comp, -3);

      if (i%100000 == 0) {
        long usedMem = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/(1024*1024);
        System.out.println("Used " + usedMem + " MB");
      }
    }
    System.out.println("-- Testing time for custom method: " + (System.currentTimeMillis() - time)/1000.0 + "s");
  }
  
  
  /**
   * Make sure that the ordering is the same for the custom binary search method
   * an a reference Collection binary sort.
   * @param iterations
   */
  private static void binarySearchTest(int iterations) {
    Random g = new Random();
    AminoAcid[] stdAA = AminoAcid.getStandardAminoAcids();
    ArrayList<Composition> keyOrder = new ArrayList<Composition>();
    
    Node node = new Node();
    long time = System.currentTimeMillis();
    for (int i=0; i<iterations; i++) {
      int comp = stdAA[(Math.abs(g.nextInt())%stdAA.length)].getComposition().getNumber();
      node.addEdge(comp, -3);
      
      Composition keyComp = new Composition(comp);
      keyOrder.add(keyComp);
    }
    
    Collections.sort(keyOrder);
    for (int i=0; i<keyOrder.size(); i++) {
      if(keyOrder.get(i).getNumber() != (int)node.getEdgeAt(i)) {
        System.out.println("FAIL at key " + node.getEdgeAt(i));
      }
    } 
    System.out.println("-- Testing time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
  }
  
  public static void main(String[] args) {
    binarySearchTest(20000);
    edgeInsertTimeMemoryTest(1000000);
  }
}
