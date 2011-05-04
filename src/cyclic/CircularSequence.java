package cyclic;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import msutil.AminoAcid;
import msutil.Composition;
import msutil.Sequence;


public class CircularSequence extends Sequence<AminoAcid>{
  private static final long serialVersionUID = 1L;

  /**
   * Constructor based on string of letters corresponding to amino acids.
   * @param sequence the sequence as a string of letters
   */
  public CircularSequence(String sequence) {
    for (int i=0; i<sequence.length(); i++) {
      AminoAcid aa = AminoAcid.getStandardAminoAcid(sequence.charAt(i));
      if (aa!=null) {
        this.add(aa);
      }
      else {
        aa = extraAATable.get(sequence.charAt(i));
        if (aa!=null) {
          this.add(aa);
        }
        else {
          System.err.println("Cannot find instantiate amino acid " + sequence.charAt(i));
          System.exit(-8);
        }
      }
    }
  }
  
  /**
   * Verify that this sequence is a subset of the peak array given as input.
   * @param peaks the peaks to check
   * @param tolerance the tolerance to incorporate the match
   * @return true if it is a subset, false otherwise.
   */
  public float getScore(ArrayList<Cluster1D> clusters, float tolerance) {
    
    float[] peaks = new float[clusters.size()];
    float[] scores = new float[clusters.size()];
    for (int i=0; i<peaks.length; i++) {
      peaks[i] = clusters.get(i).getCenter();
      scores[i] = clusters.get(i).getWeight();
    }
    
    // this sequence can be rotated up to n times
    ArrayList<Float> matches = new ArrayList<Float>();
    ArrayList<Float> bestSequence = new ArrayList<Float>(); 
    float score = 0;
    float bestScore = 0;
    for (int i=0; i<this.size(); i++) {
      
      float cumMass = 0;
      int j=0;
      matches.clear(); score = clusters.get(0).getWeight();
      for (; j<this.size()-1; j++) {
        cumMass += this.get((i+j)%this.size()).getMass();
        // search this mass in the peaks with the tolerance
        int matchIndex = find(cumMass, peaks, tolerance);
        if (matchIndex<0) {
          //System.out.println(cumMass + " not found");
          if (bestSequence.size()<matches.size()) {
            bestSequence = new ArrayList<Float>(matches);
            bestScore = score;
          }
          break;
        }
        matches.add(cumMass);
        score += scores[matchIndex];
      }
      
      // we matched everything
      if (j==this.size()-1) return score;
    }
    
    // this sequence can be rotated up to n times, and invert the sequence
    for (int i=0; i<this.size(); i++) {
      
      float cumMass = 0;
      int j=this.size()-1;
      matches.clear(); score = clusters.get(0).getWeight();
      for (; j>1; j--) {
        cumMass += this.get((i+j)%this.size()).getMass();
        // search this mass in the peaks with the tolerance
        int matchIndex = find(cumMass, peaks, tolerance); 
        if (matchIndex<0) {
          if (bestSequence.size()<matches.size()) {
            bestSequence = new ArrayList<Float>(matches); 
            bestScore = score;
          }
          break;
        }
        matches.add(cumMass);
        score += scores[matchIndex];
      }
      
      // we matched everything
      if (j==1) return score;
    }
    
    return bestScore;  
  }
  
  
  /**
   * Verify that this sequence is a subset of the peak array given as input.
   * @param peaks the peaks to check
   * @param tolerance the tolerance to incorporate the match
   * @return true if it is a subset, false otherwise.
   */
  public ArrayList<Float> isSubset(float[] peaks, float tolerance) {
    
    // this sequence can be rotated up to n times
    ArrayList<Float> matches = new ArrayList<Float>();
    ArrayList<Float> bestSequence = new ArrayList<Float>(); 
    for (int i=0; i<this.size(); i++) {
      
      float cumMass = 0;
      int j=0;
      matches.clear();
      for (; j<this.size()-1; j++) {
        cumMass += this.get((i+j)%this.size()).getMass();
        // search this mass in the peaks with the tolerance
        if (find(cumMass, peaks, tolerance)<0) {
          //System.out.println(cumMass + " not found");
          if (bestSequence.size()<matches.size()) bestSequence = new ArrayList<Float>(matches);
          break;
        }
        matches.add(cumMass);
      }
      
      // we matched everything
      if (j==this.size()-1) return matches;
    }
    
    // this sequence can be rotated up to n times, and invert the sequence
    for (int i=0; i<this.size(); i++) {
      
      float cumMass = 0;
      int j=this.size()-1;
      matches.clear();
      for (; j>1; j--) {
        cumMass += this.get((i+j)%this.size()).getMass();
        // search this mass in the peaks with the tolerance
        if (find(cumMass, peaks, tolerance)<0) {
          if (bestSequence.size()<matches.size()) bestSequence = new ArrayList<Float>(matches); 
          break;
        }
        matches.add(cumMass);
      }
      
      // we matched everything
      if (j==1) return matches;
    }
    
    return bestSequence;  
  }
  
  
  private static int find(float target, float[] items, float error) {
    int matchIndex = Arrays.binarySearch(items, target);
    if (matchIndex>=0) return matchIndex;
    
    matchIndex = -matchIndex-1;
    if (matchIndex<items.length) {
      if (items[matchIndex]-target<=error) {
        return matchIndex;
      }
      
      // check the next item if it exist
      if (matchIndex>0 && target-items[matchIndex-1]<=error) return matchIndex-1;
    }
    else {
      if (target-items[matchIndex-1]<=error) return matchIndex-1;
    }
    return -matchIndex-1;
  }
  
  
  /**
   * Verify that this sequence contains all the required breaks when aligned
   * to the given masses. The peaks created from this sequence is a superset of 
   * peaks from the list of peaks created from masses.
   * @param masses the sequence of required breaks
   * @return true if the conditions are satisfied
   */
  public boolean isSuperset(float[] masses, float tolerance) {
    return isSupersetForward(masses, tolerance) || isSupersetReverse(masses, tolerance);
  }
  
  private boolean isSupersetReverse(float[] masses, float tolerance) {
    
    // try all difference rotations
    for (int i=0; i<this.size(); i++) {
     
      int targetIndex = 0;
      float targetMass = masses[targetIndex++];
      float currentMass = 0f;
      for (int j=this.size()-1; j>=0; j--) {
        currentMass += this.get((j+i)%this.size()).getMass();
        if (currentMass < targetMass - tolerance) {
          // increment the currentMass
          continue;
        }
        
        if (currentMass <= targetMass + tolerance) {
          // increment the target mass
          if (targetIndex >= masses.length) {
            return true;
          }
          targetMass += masses[targetIndex++];
          continue;
        }

        // we missed the match
        break;
      }
    }
    
    return false;
  }
  
  private boolean isSupersetForward(float[] masses, float tolerance) {
        
    // try all difference rotations
    for (int i=0; i<this.size(); i++) {
     
      int targetIndex = 0;
      float targetMass = masses[targetIndex++];
      float currentMass = 0f;
      for (int j=0; j<this.size(); j++) {
        currentMass += this.get((j+i)%this.size()).getMass();
        if (currentMass < targetMass - tolerance) {
          // increment the currentMass
          continue;
        }
        
        if (currentMass <= targetMass + tolerance) {
          // increment the target mass
          if (targetIndex >= masses.length) {
            return true;
          }
          targetMass += masses[targetIndex++];
          continue;
        }

        // we missed the match
        break;
      }
    }
    
    return false;
  }
  
  
  public String toString() {
    StringBuffer sb = new StringBuffer();
    for (AminoAcid aa : this) {
      sb.append(aa);
      sb.append(", ");
    }
    sb.delete(sb.length()-2, sb.length());
    return sb.toString();
  }
  
  private static AminoAcid[] extraAA = 
  {
    AminoAcid.getAminoAcid('O', "Ornithine", new Composition(5, 10, 2, 1, 0)),	// modified by Sangtae
  };
  
  private static HashMap<Character,AminoAcid> extraAATable;
  static {    // initialize once
    extraAATable = new HashMap<Character,AminoAcid>();
    for (AminoAcid aa : extraAA) {
      extraAATable.put(aa.getResidue(), aa);
    }
  }

  public final static CircularSequence tyrA = new CircularSequence("VOLFPFFNQY");
  public final static CircularSequence tyrA1 = new CircularSequence("VKLFPFFNQY");
  public final static CircularSequence tyrB = new CircularSequence("VOLFPWFNQY");
  public final static CircularSequence tyrB1 = new CircularSequence("VKLFPWFNQY");
  public final static CircularSequence tyrC = new CircularSequence("VOLFPWWNQY");
  public final static CircularSequence tyrC1 = new CircularSequence("VKLFPWWNQY");
}
