package suffixarray;

import sequences.Constants;
import sequences.FastaSequence;


/**
 * This abstract class allows different formats to be searchable using a
 * SuffixArray as the database. This implementation only allows the alphabet
 * to be of sizeof(byte).
 * @author jung
 *
 */
public class SuffixArraySequence extends FastaSequence {

  
  /**
   * Constructor. The alphabet will be created dynamically according from the 
   * fasta file.
   * @param filepath the path to the fasta file.
   */
  public SuffixArraySequence(String filepath) {
    super(filepath, null);  
  }
  
  /**
   * Constructor using the specified alphabet set. If there is a letter not in
   * the alphabet.
   * @param filepath the path to the fasta file.
   * @param alphabet the specifications alphabet string. This could take the 
   *        predefined AminoAcid strings defined in this class or customized strings.
   */
  public SuffixArraySequence(String filepath, String alphabet) {
    super(filepath, alphabet,  Constants.FILE_EXTENSION);
  }
  
  /**
   * Constructor using the specified alphabet set. If there is a letter not in
   * the alphabet, it will be encoded as the TERMINATOR byte.
   * @param filepath the path to the fasta file.
   * @param alphabet the specifications alphabet string. This could take the 
   *        predefined AminoAcid strings defined in this class or customized strings.
   * @param seqExtension the extension to use for the sequence file.
   */
  public SuffixArraySequence(String filepath, String alphabet, String seqExtension) {
    super(filepath, alphabet, seqExtension);
  }
  
  /**
   * Take a ByteSequence object and make a string representation out of it.
   * @param sequence the ByteSequence object.
   * @return the translated string.
   */
  public String toString(ByteSequence sequence) {
    StringBuffer retVal = new StringBuffer(sequence.getSize());	// Switched from String to StringBuffer by sangtae
    for(int i = sequence.getSize(), index = 0; i > 0; i--, index++) {
      retVal.append(this.getCharAt(index));
    }
    return retVal.toString();
  }
  
  /**
   * This method checks whether another sequence is contained by this sequence
   * starting at a given positon.
   * @param pattern the pattern to check.
   * @param start the start position.
   * @return
   */
  public boolean contains(ByteSequence pattern, long start) {
	  if(getLCP(pattern, start) == pattern.getSize())
		  return true;
	  else
		  return false;
  }
  
  /**
   * This method returns the size of longest common prefix between pattern and a suffix of this sequence 
   * at a given positon.
   * @author sangtaekim
   * @param pattern the pattern to check.
   * @param start the start position.
   * @return
   */
  public int getLCP(ByteSequence pattern, long start) {
    long limit = Math.min(this.getSize() - start, pattern.getSize());
    int index = 0;
    for(; index < limit; index++) {
      if(pattern.getByteAt(index) != this.getByteAt(index+start)) break;
    }
    
    return index;
  }
  
  /**
   * Given a sequence translate into a byte array.
   * @param sequence the string representation.
   * @return the byte representation.
   */
  public ByteSequence toBytes(String sequence) {
    class EncodedSequence extends ByteSequence {
      private byte[] sequence;
      public EncodedSequence(byte[] sequence) {
        this.sequence = sequence;
      }
      
      public byte getByteAt(int position) {
        return this.sequence[position];
      }
      
      public int getSize() {
        return this.sequence.length;
      }
    }
    
    byte[] retSeq = new byte[sequence.length()];
    for(int i = 0; i < retSeq.length; i++) {
      if(this.isInAlphabet(sequence.charAt(i))) {
        retSeq[i] = this.toByte(sequence.charAt(i));
      }
      else {
        //retSeq[i] = sequences.Constants.TERMINATOR;
        retSeq[i] = -1;
      }
    }
    return new EncodedSequence(retSeq);
  }
  
  /**
   * Check whether this strign is fully encodable by the alphabet of this datastructure
   * @param s
   * @return
   */
  public boolean isEncodable(String s) {
    for (int index=0; index<s.length(); index++) {
      if (!isInAlphabet(s.charAt(index))) return false;
    }
    return true;
  }

}
