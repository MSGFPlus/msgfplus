package edu.ucsd.msjava.sequences;

import java.util.Collection;
import java.util.Set;

/**
 * Interface allowing access to sequence of characters. This abstract both
 * access to elements in the sequence as Characters (in original form) and Bytes
 * (in encoded form).
 * @author jung
 *
 */
public interface Sequence {

  /**
   * Return the alphabet set of this sequence as a Set of characters.
   * @return the set of characters representing the alphabet
   */
  public Collection<Character> getAlphabet();
    
  /**
   * Return the set of bytes that are valid for sequence. This is the alphabet
   * set in the form of bytes (including the terminator character, but excluding
   * un-encodable characters).
   * @return the byte alphabet set
   */
  public Set<Byte> getAlphabetAsBytes();
  
  /**
   * Returns the number of letters in the alphabet of this database.
   * @return the alphabet size including the terminator character.
   */
  public int getAlphabetSize();
 
  /**
   * Get the annotation corresponds to the given position. Annotation can be any string. 
   * For example, if the object represented is a fasta file, annotation will be the lines start with ">".
   * @param position the position to query. Annotations are mapped to certain
   *        ranges of the sequence indices, so this function will return the annotation
   *        for the subsequence that falls within this range
   * @return the annotation string corresponds to the given position. 
   */
  public String getAnnotation(long position); 
  
  /**
   * Get the encoded byte sequence at a given position.
   * No error checking for boundaries is required. 
   * @param position the position to query. 
   * @return the byte at a given position.
   */
  public byte getByteAt(long position);
 
  /**
   * Get a slice of the sequence as a byte array representation.
   * @param start the start index
   * @param end the end index (exclusive)
   * @return the byte array of the slice
   */
  public byte[] getBytes(int start, int end);
  
  /**
   * Retrieve the original character in the sequence at the given position.
   * No error checking for boundaries is required.
   * @param position the location of the interested character.
   * @return the character at the given position.
   */
  public char getCharAt(long position);
  
  /**
   * Get the unique identified for this sequence.
   * @return the unique identifier used for the suffix array to verified that the tree was built on the same sequence.
   */
  public int getId();
  
  /** 
   * Get the complete entry that corresponds to the given position. 
   * A chunk is a contiguous sequence that shares an annotation (e.g. Protein).
   * @param position the position to query.
   * @return the entry as a string corresponds to the given position. 
   */
  public abstract String getMatchingEntry(long position);
  
  /** 
   * Get the complete entry that corresponds to the given fasta entry header. 
   * A chunk is a contiguous sequence that shares an annotation (e.g. Protein).
   * @param name the header of the fasta entrry.
   * @return the entry as a string. 
   */
  public abstract String getMatchingEntry(String name);
  
  /**
   * The size of this sequence. All indexes [0, getSize()) should have valid characters.
   * @return the size of this sequence.
   */
  public long getSize();
  
  /**
   * Check whether the given character is part of the (specified) alphabet.
   * @param c the character to check
   * @return the membership of the char in the alphabet 
   */
  public boolean isInAlphabet(char c);
  
  /**
   * A quick way to find out whether the given position corresponds to terminating
   * character
   * @param position the position to inquire about
   * @return true if the given position is a terminator, false otherwise
   */
  public boolean isTerminator(long position);

  /**
   * Check that the character at this position is not a terminator and it is in
   * the alphabet.
   * @param position the position
   * @return the truth value of the statement above.
   */
  public boolean isValid(long position);
  
  /**
   * Translates the given character to the binary representation.
   * @param c the character to convert.
   * @return the binary representation if available, but this might fail if the
   *         character is not in the alphabet.
   */
  public byte toByte(char c);
  
  /**
   * Take a byte and reverse translate it to the original string representation.
   * Note that some bytes might represent more than one character. An arbitrary
   * character is returned. To find out what the original char was, the getCharAt
   * method should be called instead
   * @param b the byte.
   * @return the String representation of the given byte.
   */
  public char toChar(byte b);
  
  /**
   * Translates from a byte sequence to a character sequence.
   * @param sequence the array of bytes to translate.
   * @return the string representation of the given sequence.
   */
  public String toString(byte[] sequence);

  /**
   * Returns the starting position of the protein covered by the coordinate
   * given by the parameter
   * @param position any location of the protein
   * @return the starting position
   */
  public long getStartPosition(long position);
  
  /**
   * Get a slice of the string by the given coordinates.
   * @param start the starting position (inclusive)
   * @param end the ending position (exclusive)
   * @return the string representation of the given subsequence.
   */
  public String getSubsequence(long start, long end);
  
  
}
