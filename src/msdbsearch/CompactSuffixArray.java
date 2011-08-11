package msdbsearch;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.IntBuffer;
import java.util.Arrays;

import sequences.Constants;
import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;
import suffixarray.SuffixFactory;

public class CompactSuffixArray extends SuffixArray {

	public CompactSuffixArray(SuffixArraySequence sequence) {
		super(sequence);
	}

	/**
	 * Helper method that creates the suffixFile.
	 * @param sequence the Adapter object that represents the database (text).
	 * @param suffixFile the output file.
	 */
	@Override
	protected void createSuffixArrayFile(SuffixArraySequence sequence, String suffixFile) {
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
					// JAVA 1.5 code
					int[] tempArray = new int[this.size+INCREMENT_SIZE];
					for(int i = 0; i < size; i++)     tempArray[i] = this.items[i];
					this.items = tempArray;
				}
				/* JAVA 1.6 code
          		this.items = Arrays.copyOf(this.items, this.size+INCREMENT_SIZE);
				 */
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
		{
			denominator *= hashBase;
		}

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
		for(int i=BUCKET_SIZE-1, j=0, limit = (int)sequence.getSize(); j < limit; i++, j++) {

			// print progress
			if(j % 1000001 == 0)  System.out.printf("Suffix creation: %.2f%% complete.\n", j*100.0/sequence.getSize());

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
			DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(suffixFile)));
			out.writeInt((int)sequence.getSize());
			out.writeInt(sequence.getId());      
			SuffixFactory.Suffix prevBucketSuffix = null;
			byte[] neighboringLcps = new byte[(int)sequence.getSize()];         // the computed neighboring lcps
			int order = 0;
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
					out.writeInt(first.getIndex());
					neighboringLcps[order++] = lcp; 
					SuffixFactory.Suffix prevSuffix = first;

					for(int j = 1; j < sortedSuffixes.length; j++) {
						SuffixFactory.Suffix thisSuffix = sortedSuffixes[j];
						//store the information
						out.writeInt(thisSuffix.getIndex());
						neighboringLcps[order++] = thisSuffix.getLCP(prevSuffix, BUCKET_SIZE);
						prevSuffix = thisSuffix;
					}
					prevBucketSuffix = sortedSuffixes[0];
				}
			}

			// compute the leftMiddle and middleRight lcps
			System.out.println("Computing the parameterized lcp arrays..");
			out.write(neighboringLcps);	// Sangtae
			out.flush(); 
			out.close();
		}
		catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return;
	}	
	@Override
	protected int readSuffixArrayFile(String suffixFile) {
		//		System.out.println("SAForMSGFDB Reading " + suffixFile);
		try {
			// read the first integer which encodes for the size of the file
			DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(suffixFile)));
			size = in.readInt();
			// the second integer is the id
			int id = in.readInt();

			int[] indexArr = new int[size];
			for(int i=0; i<indexArr.length; i++)
				indexArr[i] = in.readInt();
			indices = IntBuffer.wrap(indexArr).asReadOnlyBuffer();

			int sizeOfLcps = size;
			// skip leftMiddleLcps and middleRightLcps
			long totalBytesSkipped = 0;
			while(totalBytesSkipped < 2*sizeOfLcps)
			{
				long bytesSkipped = in.skip(2*sizeOfLcps-totalBytesSkipped);
				if(bytesSkipped == 0)
				{
					System.out.println("Error while reading suffix array: " + totalBytesSkipped + "!=" +2*sizeOfLcps);
					System.exit(-1);
				}
				totalBytesSkipped += bytesSkipped;
			}
			if(totalBytesSkipped != 2*sizeOfLcps)
			{
				System.out.println("Error while reading suffix array: " + totalBytesSkipped + "!=" +2*sizeOfLcps);
				System.exit(-1);
			}
			// read neighboringLcps
			byte[] neighboringLcpArr = new byte[sizeOfLcps];
			in.read(neighboringLcpArr);
			neighboringLcps = ByteBuffer.wrap(neighboringLcpArr).asReadOnlyBuffer();
			in.close();

			return id;
		}
		catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

		return 0;
	}		
}
