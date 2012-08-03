package edu.ucsd.msjava.suffixarray;


/**
 * This abstract class allows the query of the suffix array.
 * @author jung
 *
 */
public abstract class ByteSequence implements Comparable<ByteSequence> {

	public static final int MAX_COMPARISON_LENGTH = Byte.MAX_VALUE;
	// maximum number of characters to print for this sequence
	private final int PRINT_LIMIT = 80;

	/**
	 * Get the byte value at the given index.  
	 * @param index the index to retrieve the value from.
	 * @return the byte at position index.
	 */
	public abstract byte getByteAt(int index); 


	/**
	 * Get the size of this sequence.
	 * @return the size of this sequence.
	 */
	public abstract int getSize();


	/**
	 * Lexographic compare that break once it finds the tie breaker index before
	 * hashSize.
	 * @param other other suffix to compare against.
	 */
	public int compareTo(ByteSequence other) {
		return compareTo(other, 0);
	}


	/**
	 * Returns a copy of the byte sequence encoded by this array.
	 * @return a byte array.
	 */
	public byte[] getSequence() {
		byte[] sequence = new byte[getSize()];
		for(int i=0, limit=getSize(); i < limit; i++)     sequence[i] = getByteAt(i);
		return sequence;
	}


	/**
	 * Lexographic compare that break once it finds the tie breaker index before
	 * hashSize.
	 * @param other other suffix to compare against.
	 * @param start start comparing from this offset.
	 * @return a positive number if this > other, a negative if other > this and 0
	 *         if they are equal. The longest common prefix length can be retrieved
	 *         by taking absolute value of the return value minus 1.
	 */
	public int compareTo(ByteSequence other, int start) {
		int limit = Math.min(this.getSize(), other.getSize());
		if(limit > MAX_COMPARISON_LENGTH)
			limit = MAX_COMPARISON_LENGTH;
		int offset = start;
		for(; offset < limit; offset++) {
			byte thisByte = this.getByteAt(offset);
			byte otherByte = other.getByteAt(offset);
			if(thisByte > otherByte)       return offset+1;
			else if(otherByte > thisByte)  return -offset-1;
		}
		// the longer one is the greater one
		if(this.getSize() > other.getSize())    return offset+1;
		if(other.getSize() > this.getSize())    return -offset-1;
		return 0;
	}


	/**
	 * Calculates the index of the longest common prefix of two given suffixes.
	 * @param other the suffix to compare against.
	 * @param start the starting position to start comparing.
	 * @return the number that is returned is the number of positions in which
	 *         the prefixes of the two objects are equal. 0 means that nothing
	 *         is common.
	 */
	public byte getLCP(ByteSequence other, int start) {
		int limit = Math.min(this.getSize(), other.getSize());
		if(limit > Byte.MAX_VALUE)
			limit = Byte.MAX_VALUE;
		int offset = start;
		for(; offset < limit; offset++) {
			byte thisByte = this.getByteAt(offset);
			byte otherByte = other.getByteAt(offset);
			if(thisByte > otherByte || otherByte > thisByte)   
				return (byte)offset;
		}
		return (byte)offset;
	}

	/**
	 * Calculates the index of the longest common prefix of two given suffixes.
	 * @param other the suffix to compare against.
	 * @param start the starting position to start comparing.
	 * @return the number that is returned is the number of positions in which
	 *         the prefixes of the two objects are equal. 0 means that nothing
	 *         is common.
	 */
	public int getLCPInt(ByteSequence other, int start) {
		int limit = Math.min(this.getSize(), other.getSize());
		int offset = start;
		for(; offset < limit; offset++) {
			byte thisByte = this.getByteAt(offset);
			byte otherByte = other.getByteAt(offset);
			if(thisByte > otherByte || otherByte > thisByte)   
				return offset;
		}
		return offset;
	}


	/**
	 * Overloaded method to get the LCP.
	 * @param other the other suffix to compare.
	 * @return the lcp of this and other.
	 */
	public byte getLCP(ByteSequence other) {
		return getLCP(other, 0);
	}

	/**
	 * Overloaded method to get the LCP.
	 * @param other the other suffix to compare.
	 * @return the lcp of this and other.
	 */
	public int getLCPInt(ByteSequence other) {
		return getLCPInt(other, 0);
	}


	/**
	 * Return the string representation of this sequence.
	 * @return the string representing the string of this sequence.
	 */
	public String toString() {
		String retVal = "";
		for(int i = 0, limit = Math.min(getSize(), PRINT_LIMIT); i < limit; i++) {
			retVal += getByteAt(i) + " ";
		}
		return retVal;
	}

}
