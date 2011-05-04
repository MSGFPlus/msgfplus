/**
 * 
 */
package mstag;

import java.util.*;

import parser.MgfSpectrumParser;
import msutil.*;



/**
 * Contains the static methods to generate tags given a spectrum.
 * @author jung
 *
 */
public class Tagger {
  
  private final static AminoAcid[] stdAminoAcids;
  private final static float maxMass;
  static {
    // ignore I, and use L instead
    stdAminoAcids = AminoAcid.getStandardAminoAcids();
    Arrays.sort(stdAminoAcids);
    maxMass = stdAminoAcids[stdAminoAcids.length-1].getMass();
  }
  
  
  /**
   * A wrapper of the binary search algorithm that accounts for tolerance and
   * returns all matches.
   * @param mass the target mass.
   * @param tolerance how much slack to allow for the matches.
   * @return the array of amino acids 
   */
  private final static AminoAcid[] getAminoAcids(Float mass, float tolerance) {
    AminoAcid target = AminoAcid.getCustomAminoAcid('X', mass);
    int index = Arrays.binarySearch(stdAminoAcids, target);
    if(index < 0)  index = -index - 1;
    
    LinkedList<AminoAcid> results = new LinkedList<AminoAcid>();

    // need to hover in the neighborhood to truly see whether there are no matches
    while(index < stdAminoAcids.length && stdAminoAcids[index].getMass() <= mass + tolerance) {
      index++;
    }
    index--;
    
    while(index >= 0 && stdAminoAcids[index].getMass() >= mass - tolerance) {
      results.addFirst(stdAminoAcids[index]);
      index--;
    }
    return results.toArray(new AminoAcid[0]);
  }
  
  
  /**
   * This class only has static methods.
   */
  private Tagger() {
  }

  
  /**
   * Main method. Testing code.
   * @param args the command line paremeters.
   * @throws Exception 
   */
  public static void main(String[] args) throws Exception{

    float cutOffFraction = 0.5f;
    float ppmTolerance = 10f;
    
    String spectrumFile = "/Users/jung/Research/USTags/centroidedMGF/spec1845.mgf";
    
    SpectraContainer sc = new SpectraContainer(spectrumFile, new MgfSpectrumParser());
    Reshape wf = new RankFilter(380);
    Spectrum s = wf.apply(sc.get(0));
    
    System.out.println("Parent mass " + s.getParentMass());
    
    // log the intensities of the peaks
    for(Peak p : s) {
      p.setIntensity((float)Math.log(p.getIntensity()));
    }
    
    System.out.println("Peak count before fitlering "+s.size());
    ArrayList<Tag> results = makeTags(s, ppmTolerance, cutOffFraction);
    System.out.println(results.size());
    
    // print top 10 results
    int rank = 1;
    for(Tag tag : results) {
      System.out.println(rank + "\t" + tag);
      rank += 1;
      if(rank > 10)     break;
    }
    
    rank = 1;
    for(Tag tag : results) {
      String seqStr = tag.getTagStr();
      if(seqStr.equals("GFTFSFPAS")) {
        System.out.println(rank + "\t" + tag);
      }
      rank += 1;
    }
    
  }

  
 
  public static ArrayList<Tag> makeTags(Spectrum s, float ppmTolerance, float cutOffFraction) {
  
    // this keeps track of the peaks we need to examine
    LinkedList<Peak> peakBuffer = new LinkedList<Peak>();
    
    // this is the dp table
    Cell[] table = new Cell[s.size()];
    for(int i = table.length-1; i >= 0; i--) {
      table[i] = new Cell(s.get(i));
      table[i].setScore(s.get(i).getIntensity());     // initialize
    }
  
    float maxScore = Float.NEGATIVE_INFINITY;
    
    // fill the table
    int peakIndex = 0;
    for(Peak current : s) {
      // remove unreachable peaks from the buffer
      while(peakBuffer.size() > 0) {
        if(peakBuffer.getFirst().getMass() < current.getMass() - maxMass) {
          peakBuffer.removeFirst();
        }
        else {
          break;
        }
      }
      
      for(Peak previous : peakBuffer) {
        float massDiff = current.getMass() - previous.getMass();
        AminoAcid[] matches = getAminoAcids(massDiff, current.toUnitTolerance(ppmTolerance));
        for(AminoAcid aa : matches) {
          Cell previousCell = table[previous.getIndex()];
          float thisDelta = current.getIntensity();   // this is the actual score 
          float thisScore = previousCell.getScore() + thisDelta;
          table[peakIndex].addLink(thisScore, previousCell, aa);
          if(thisScore > maxScore)     maxScore = thisScore;
        }
        
      }

      current.setIndex(peakIndex);
      peakBuffer.addLast(current);
      peakIndex++;
    }
   
    // retrieve the sequences recursively
    float cutOff = maxScore * cutOffFraction; 
    ArrayList<Tag> results = new ArrayList<Tag>();
    for(Cell c: table) {
      if(c.score >= cutOff) {
        Stack<Peak> peaks = new Stack<Peak>();
        peaks.push(c.peak);
        retrieve(c, new Stack<AminoAcid>(), peaks, c.peak.getIntensity(), cutOff, s.getParentMass() - c.peak.getMass(), results);
      }
    }
    
    Collections.sort(results);  
    return results; 
  }
  
  
  
  /**
   * Retrieve the sequences that satisfy a given cutoff.
   * @param target
   * @param sequence
   * @param peaks
   * @param cumScore
   * @param cutOff
   * @param rightMass
   * @param results
   */
  private static void retrieve(Cell target,
                               Stack<AminoAcid> sequence, 
                               Stack<Peak> peaks, 
                               float cumScore, 
                               float cutOff,
                               float rightMass,
                               ArrayList<Tag> results) {
        
    // if the sequence has enough score add it
    float thisScore = cumScore + target.peak.getIntensity(); 
    if(thisScore >= cutOff) {
      Peptide tag = new Peptide(sequence.toArray(new AminoAcid[sequence.size()]));
      Peak[] finalPeaks = peaks.toArray(new Peak[peaks.size()]);
      Mass left = new Mass(peaks.get(0).getMass());
      Mass right = new Mass(rightMass);
      results.add(new Tag(left, tag, right, finalPeaks, thisScore));
    }
    
    for(Cell.Link link : target.getLinks()) {  
  
      Cell previous = link.previous;      
      
      // maximum achievable score
      float bestScore = previous.score + cumScore;
      if(bestScore >= cutOff) {

        // prepare to recurse
        float newCumScore = target.peak.getIntensity() + cumScore;
        sequence.push(link.aa);
        peaks.push(target.peak);
        
        // recurse
        retrieve(previous,
                 sequence,
                 peaks,
                 newCumScore,
                 cutOff,
                 rightMass,
                 results);
        
        // restore to pre recursion state
        sequence.pop();
        peaks.pop();
      }  
    }
    
  }
  
  
  
  
  /**
   * The cell that holds in information for the dynamic programming table.
   * @author jung
   *
   */
  private static class Cell {
    
    private class Link implements Comparable<Link>{
      private float score;   // score ending at this link 
      private Cell previous; // how to backtrack
      private AminoAcid aa;  // what amino acid corresponds this
      
      private Link(float score, Cell previous, AminoAcid aa) {
        this.score = score;
        this.previous = previous;
        this.aa = aa;
      }
      
      /**
       * The score gives the order for a links. Use reverse order.
       * @param other the other object to compare to.
       * @return 1 if other > this, -1 if this > other, 0 otherwise.
       */
      public int compareTo(Link other) {
        if(this.score > other.score)   return -1;
        if(other.score > this.score)   return 1;
        return 0;
      }
      
      
    }
    
    // the best score of this cell
    private float score;
    
    // each cell corresponds to a peak
    private Peak peak;
    
    // backtracking options
    private PriorityQueue<Link> links;
    
    private Cell(Peak peak) {
      score = Float.NEGATIVE_INFINITY;
      links = new PriorityQueue<Link>();
      this.peak = peak;
    }
    
    public void setScore(float score) {
      this.score = score;
    }
    
    public float getScore() {
      return score;
    }
    
    public Collection<Link> getLinks() {
      return links;
    }
    
    private void addLink(float score, Cell previous, AminoAcid aa) {
      if(this.score < score)      this.score = score;
      links.add(new Link(score, previous, aa));
    }
    
    
  }
}
