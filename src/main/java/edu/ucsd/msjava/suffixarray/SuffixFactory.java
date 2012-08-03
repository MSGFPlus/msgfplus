package edu.ucsd.msjava.suffixarray;

import edu.ucsd.msjava.sequences.Sequence;



/**
 * SuffixFactory and Suffix classes. This class will allow the creation of 
 * light weight suffix objects given a long sequence in the form of an Adapter
 * object.
 * @author jung
 *
 */
public class SuffixFactory {


	/**
	 * Class that represents a Suffix object.
	 * @author jung
	 *
	 */
	public class Suffix extends ByteSequence {

		// the index of this suffix
		private int index;
		// modified by Sangtae to save memory
		//    private int size;


		/**
		 * Constructor.
		 * @param index the starting index of the suffix.
		 */
		public Suffix(int index) {
			this.index = index;
			//      this.size = (int)sequence.getSize() - index;
		}


		public int getSize() {
			//      return this.size;
			return (int)sequence.getSize() - index;
		}


		public byte getByteAt(int index) {
			return sequence.getByteAt(this.index+index);
		}


		/**
		 * Getter method.
		 * @return the index of this suffix.
		 */
		public int getIndex() {
			return this.index;
		}
	}


	// modified by Sangtae
	// holds the sequences
	//	private SuffixArraySequence sequence;
	private Sequence sequence;

	/**
	 * Constructor.
	 * @param sequence the sequence object to create the suffixes from.
	 */
	public SuffixFactory(Sequence sequence) {
		this.sequence = sequence;
	}


	/**
	 * Factory method that creates a new suffix object from a sequence.
	 * @param index the starting index of the suffix.
	 * @return the suffix object
	 */
	public Suffix makeSuffix(int index) {
		return new Suffix(index);
	}


	/**
	 * Get the longest common prefix count for 2 given suffixes. 
	 * @param o1 one of the objects.
	 * @param o2 the other object.
	 * @param offset the number of indexes to skip when calculating the LCP.
	 * @return the number of positions in which these 2 suffixes are in common
	 *         or offset (if this number is greater). 
	 */
	public int getLCP(Suffix o1, Suffix o2, int offset) {
		return o1.getLCP(o2, offset);
	}


	/**
	 * Overloaded method where the offset is 0.
	 * @param o1 one of the objects.
	 * @param o2 the other object.
	 * @return refer to the documentation of the other method.
	 */
	public int getLCP(Suffix o1, Suffix o2) {
		return o1.getLCP(o2, 0);
	}

}
