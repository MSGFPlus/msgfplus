package msdbsearch;

import java.io.*;
import java.util.*;
import sequences.Constants;
import suffixarray.SuffixFactory;

/**
 * SuffixArray class for fast exact matching.
 * @author Sangtae Kim
 *
 */
public class CompactSuffixArray {

	/***** CONSTANTS *****/
	// The default extension of a suffix array file.
	protected static final String EXTENSION_INDICES = ".sarr";
	protected static final String EXTENSION_NLCPS = ".snlcp";

	// the size of the bucket for the suffix array creation
	protected static final int BUCKET_SIZE = 5;

	// the size of an int primitive type in bytes
	protected static final int INT_BYTE_SIZE = Integer.SIZE / Byte.SIZE;

	/***** MEMBERS *****/
	// the indices of the sorted suffixes
//	DataInputStream indices;
	private final File indexFile;
	
	// precomputed LCPs of neighboring suffixes
//	DataInputStream neighboringLcps;
	private final File nlcpFile;
	
	// the sequence representing all the suffixes
	private CompactFastaSequence sequence;

	// the class that generates suffixes from the given adapter
	private SuffixFactory factory;

	// the number of suffixes in this array
	private int size;

	// the number of distinct peptides
	private final int maxPeptideLength;
	private int[] numDisinctPeptides;
	
	/**
	 * Constructor that attempts to read the suffix array from the provided file.
	 * @param sequence the sequence object.
	 */
	public CompactSuffixArray(CompactFastaSequence sequence, int maxPeptideLength) {
		// infer the suffix array file from the sequence.
		this.sequence = sequence;
		this.factory = new SuffixFactory(sequence);
		indexFile = new File(sequence.getBaseFilepath() + EXTENSION_INDICES);
		nlcpFile = new File(sequence.getBaseFilepath() + EXTENSION_NLCPS);

		// create the file if it doesn't exist.
		if(!indexFile.exists() || !nlcpFile.exists()) {           
			createSuffixArrayFiles(sequence, indexFile, nlcpFile);
		}

		// check the ids of indexFile and nlcpFile
		int id = checkID();

		// check that the files are consistent
		if(id != sequence.getId()) {
			System.err.println("Suffix array files are not consistent: " + indexFile + ", " + nlcpFile);
			System.err.println("Please recreate the suffix array file.");
			System.exit(-1);
		}
		
		this.maxPeptideLength = maxPeptideLength;
		computeNumDistinctPeptides();
	}

	public File getIndexFile()	{	return this.indexFile; }
	public File getNeighboringLcpFile()	{	return this.nlcpFile; }
	public CompactFastaSequence getSequence()	{ return sequence; }
	
	public int getSize()
	{
		return size;
	}

	public int getNumDistinctPeptides(int length)	
	{
		// no boundary check
		return numDisinctPeptides[length];
	}
	
	public String getAnnotation(long index)
	{
		return sequence.getAnnotation(index);
	}	
	
	private void computeNumDistinctPeptides()
	{
		numDisinctPeptides = new int[maxPeptideLength+2];
		try {
			DataInputStream neighboringLcps = new DataInputStream(new BufferedInputStream(new FileInputStream(nlcpFile)));
			int size = neighboringLcps.readInt();
			neighboringLcps.readInt();	// skip id
			
			for(int i=0; i<size; i++)
			{
				byte lcp = neighboringLcps.readByte();
				for(int l=lcp+1; l<numDisinctPeptides.length; l++)
				{
					numDisinctPeptides[l]++;
				}
			}
			neighboringLcps.close();
		}
		catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	/**
	 * Helper method that initializes the suffixArray object from the file.
	 * Initializes indices, leftMiddleLcps, middleRightLcps and neighboringLcps.
	 * @param suffixFile the suffix array file.
	 * @return returns the id of this file for consistency check.
	 */
	private int checkID() {
		//		System.out.println("SAForMSGFDB Reading " + suffixFile);
		try {
			DataInputStream indices = new DataInputStream(new BufferedInputStream(new FileInputStream(indexFile)));
			// read the first integer which encodes for the size of the file
			int sizeIndexFile = indices.readInt();
			// the second integer is the id
			int idIndexFile = indices.readInt();

			DataInputStream neighboringLcps = new DataInputStream(new BufferedInputStream(new FileInputStream(nlcpFile)));
			int sizeNLcp = neighboringLcps.readInt();
			int idNLcp = neighboringLcps.readInt();
			
			indices.close();
			neighboringLcps.close();
			
			if(sizeIndexFile == sizeNLcp && idIndexFile == idNLcp)
				return idIndexFile;
		}
		catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

		return 0;
	}		
	
	/**
	 * Helper method that creates the suffixFile.
	 * @param sequence the Adapter object that represents the database (text).
	 * @param suffixFile the output file.
	 */
	private void createSuffixArrayFiles(CompactFastaSequence sequence, File indexFile, File nlcpFile) {
		System.out.println("Creating the suffix array indexed file... Size: " + sequence.getSize());

		// helper local class
		class Bucket {
			// how much to increment once we reach the maximum occupancy for a bucket
			private static final int INCREMENT_SIZE = 10;
			private int[] items;
			private int size;

			/**
			 * Constructor.
			 */
			public Bucket() {
				this.items = new int[10];
				this.size = 0;
			}

			/**
			 * Add item to the bucket.
			 * @param item the item to add.
			 */
			public void add(int item) {
				if(this.size >= items.length) {
					this.items = Arrays.copyOf(this.items, this.size+INCREMENT_SIZE);
				}
				this.items[this.size++] = item;
			} 

			/**
			 * Get a sorted version of this bucket.
			 * @return
			 */
			public SuffixFactory.Suffix[] getSortedSuffixes() {
				SuffixFactory.Suffix[] sa = new SuffixFactory.Suffix[this.size];
				for(int i = 0; i < this.size; i++) {
					sa[i] = factory.makeSuffix(this.items[i]);
				}
				Arrays.sort(sa);
				return sa;
			}
		}

		// the size of the alphabet to make the hashes
		int hashBase = sequence.getAlphabetSize();
		if(hashBase > 30)
		{
			System.err.println("Suffix array construction failure: alphabet size is too large: " + sequence.getAlphabetSize());
			System.exit(-1);
		}

		// this number is to efficiently calculate the next hash
		int denominator = 1;
		for(int i=0; i < BUCKET_SIZE-1; i++)
			denominator *= hashBase;

		// the number of buckets  required to encode for all hashes
		int numBuckets = denominator*hashBase;

		// initial value of the hash
		int currentHash = 0;
		for(int i=0; i < BUCKET_SIZE-1; i++) {    
			currentHash = currentHash*hashBase + sequence.getByteAt(i);
		}

		// the main array that stores the sorted buckets of suffixes
		Bucket[] bucketSuffixes = new Bucket[numBuckets];

		// main loop for putting suffixes into the buckets
		for(int i=BUCKET_SIZE-1, j=0; j < (int)sequence.getSize(); i++, j++) {
			// print progress
			if(j % 1000001 == 0)  System.out.printf("Suffix creation: %.2f%% complete.\n", j*100.0/sequence.getSize());

//			System.out.println(j);
			// quick wait to derive the next hash, since we are reading the sequence in order 
			byte b = Constants.TERMINATOR;
			if (i<sequence.getSize()) b = sequence.getByteAt(i);
			currentHash = (currentHash%denominator)*hashBase+b;

			// first bucket at this position
			if(bucketSuffixes[currentHash] == null)    bucketSuffixes[currentHash] = new Bucket();

			// insert suffix
			bucketSuffixes[currentHash].add(j);
		}

		try {
			DataOutputStream indexOut = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(indexFile)));
			DataOutputStream nLcpOut = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(nlcpFile)));
			indexOut.writeInt((int)sequence.getSize());
			indexOut.writeInt(sequence.getId());      
			nLcpOut.writeInt((int)sequence.getSize());
			nLcpOut.writeInt(sequence.getId());
			SuffixFactory.Suffix prevBucketSuffix = null;
//			byte[] neighboringLcps = new byte[(int)sequence.getSize()];         // the computed neighboring lcps
			for(int i=0; i < bucketSuffixes.length; i++) {

				// print out progress
				if(i % 100000 == 99999)   System.out.printf("Sorting %.2f%% complete.\n", i*100.0/bucketSuffixes.length);

				if(bucketSuffixes[i] != null) {

					SuffixFactory.Suffix[] sortedSuffixes = bucketSuffixes[i].getSortedSuffixes();

					SuffixFactory.Suffix first = sortedSuffixes[0];
					byte lcp = 0;
					if(prevBucketSuffix != null) {
						lcp = first.getLCP(prevBucketSuffix);
					}
					// write information to file
					indexOut.writeInt(first.getIndex());
//					neighboringLcps[order++] = lcp;
					nLcpOut.writeByte(lcp);
					SuffixFactory.Suffix prevSuffix = first;

					for(int j = 1; j < sortedSuffixes.length; j++) {
						SuffixFactory.Suffix thisSuffix = sortedSuffixes[j];
						//store the information
						indexOut.writeInt(thisSuffix.getIndex());
//						neighboringLcps[order++] = thisSuffix.getLCP(prevSuffix, BUCKET_SIZE);
						nLcpOut.writeByte(thisSuffix.getLCP(prevSuffix, BUCKET_SIZE));
						prevSuffix = thisSuffix;
					}
					prevBucketSuffix = sortedSuffixes[0];
					bucketSuffixes[i] = null;	// deallocate the memory
				}
			}

			bucketSuffixes = null;
			
			indexOut.flush();
			indexOut.close();
			
			nLcpOut.flush();
			nLcpOut.close();
			
			// Do not compute Llcps and Rlcps
		}
		catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return;
	}

	@Override
	public String toString() {
		String retVal = "Size of the suffix array: " + this.size + "\n";
//		int rank = 0;
//		while(indices.hasRemaining()) {
//			int index = indices.get();
//			int lcp = this.neighboringLcps.get(rank);
//			retVal += rank + "\t" + index + "\t" + lcp + "\t" + sequence.toString(factory.makeSuffix(index).getSequence()) + "\n";
//			rank++;
//		}
//		indices.rewind();        // reset marks after iteration
//		neighboringLcps.rewind();
		return retVal;
	}
}
