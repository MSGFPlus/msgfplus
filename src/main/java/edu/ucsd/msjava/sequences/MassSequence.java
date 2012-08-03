package edu.ucsd.msjava.sequences;

public interface MassSequence extends Sequence {

  /**
   * This is a special method to handle protein fasta sequences in which allows
   * the query of a mass of an amino at a certain index.
   * @param index the index of the item in which mass we want to know
   * @return the integer mass of the amino acid at the given position, 0 if the
   *         position corresponds to a TERMINATOR or unknown amino acid.
   */
  public int getIntegerMass(long index);
  
  /**
   * Calculates the mass of a segment of this sequence.
   * @param start the start of the segment (inclusive)
   * @param end the end of segment (exclusive)
   * @return the integer mass of the given segment. If there are unknown amino
   * acids in the segment, their masses will be treated as 0. 
   */
  public int getIntegerMass(long start, long end);
 
  /**
   * Checks whether this position can be translated into a mass.
   * @param position the position to check
   * @return true if this has a mass, false otherwise.
   */
  public boolean hasMass(long position);
  
}
