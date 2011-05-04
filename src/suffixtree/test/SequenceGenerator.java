package suffixtree.test;

import java.util.ArrayList;
import java.util.Random;

import sequences.MassSequence;
import sequences.Sequence;
import suffixtree.Constants;

public class SequenceGenerator {
  

  /**
   * Quickly return the position of the next invalid or terminator character.  
   * @param start the position to start looking from
   * @return the position of the terminating or invalid character
   */
  private static long getNextEnd(Sequence s, long start) {
    for (long i=start+1; i<s.getSize(); i++) {
      if (s.isTerminator(i)) return i;
    }
    return s.getSize();
  }
  
  
  /**
   * Generate a mass set representing a substring of this sequence starting from
   * the specified position and satisfying the other parameters.
   * @param s the sequence to use
   * @param start the start index where to generate the sequence
   * @param lower the minimum mass of the returned sequence
   * @param upper the maximum mass of the returned sequence
   * @return the array list of integer masses or null if such sequence cannot be
   *         constructed from start.
   */
  public static ArrayList<Integer> generateIntegerMassSet(MassSequence s, long start, int lower, int upper) {
    if (s.isTerminator(start)) return null;
    
    long end = getNextEnd(s, start);
    
    if (end > start) {
      ArrayList<Integer> results =  new ArrayList<Integer>();
      int cumMass = 0;
      int mass;
      while (start<end) {
        mass = s.getIntegerMass(start++);
        cumMass += mass;
        if (cumMass > upper) break;
        results.add(mass);
      }
      if (lower <= cumMass) return results;
    }
    return null;
  }
  
  
  /**
   * Generate a gapped mass set of the amino acid sequence starting
   * from a given location with the specified parameters.
   * @param s the sequence to use
   * @param start the start postions
   * @param lower the minimum mass of the returned sequence
   * @param upper the maximum mass of the returned sequence
   * @param maxGapMass the largest gap possible as Daltons
   * @param r the random number generator object
   * @return the arraylist of integer masses if there exist a sequence that
   *         statisfy the parameters, or null otherwise.
   */
  public static ArrayList<Integer> generateRandomIntegerMassSet(MassSequence s, long start, int lower, int upper, int maxGapMass, Random r) {
    if (s.isTerminator(start)) return null;
    
    long end = getNextEnd(s, start);
    
    if (end > start) {
      int superTotal = 0;
      ArrayList<Integer> results = new ArrayList<Integer>();
      while (start<end) {
        int topMass = r.nextInt(maxGapMass); 
        int cumMass = s.getIntegerMass(start++);
        
        while (start<end) {
          if (cumMass+s.getIntegerMass(start)<=topMass) {
            cumMass += s.getIntegerMass(start++);
          }
          else break;
        }
        superTotal += cumMass;
        if (superTotal>upper) break;
        results.add(cumMass);
      }
      if (lower <= superTotal) return results;
    }
    return null;
  }
  
  
  /**
   * Generate a given number of random correct queries.
   * @param s the sequence to use
   * @param queries the data structure to store the results
   * @param count the number of queries to generate
   */
  public static void generateRandomCorrectQueries(MassSequence s, ArrayList<ArrayList<Integer>> queries, int count) {
    Random r = new Random();
    
    for (int i=0; i<count; i++) {
      int start = r.nextInt((int)s.getSize());
      ArrayList<Integer> gapQuery = generateRandomIntegerMassSet(s, start, Constants.MIN_QUERY_MASS, Constants.MAX_QUERY_MASS, Constants.MAX_GAP_MASS, r);   
      if (gapQuery==null) {
        i--;
      }
      else {
        queries.add(gapQuery);
      }
    }
  }
  
  /**
   * Generate all correct queries for all starting positions.
   * @param s the sequence to use
   * @param queries the data structure to store the results
   * @param count the number of queries to generate
   */
  public static void generateAllCorrectQueries(MassSequence s, ArrayList<ArrayList<Integer>> queries) {
    Random r = new Random();
    
    for (int start=0; start<(int)s.getSize(); start++) {
      ArrayList<Integer> gapQuery = generateRandomIntegerMassSet(s, start, Constants.MIN_QUERY_MASS, Constants.MAX_QUERY_MASS, Constants.MAX_GAP_MASS, r);   
      if (gapQuery!=null) {
        queries.add(gapQuery);
      }
    }
  }
}
