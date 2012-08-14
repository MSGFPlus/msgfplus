package edu.ucsd.msjava.suffixarray;

import java.io.*;
import java.nio.*;
import java.nio.channels.FileChannel;
import java.util.*;

import edu.ucsd.msjava.msdbsearch.CompactSuffixArray;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.sequences.Constants;


/**
 * SuffixArray class for fast exact matching.
 * @author Sangtae Kim
 *
 */
public class SuffixArray {



	/***** CONSTANTS *****/
	// The default extension of a suffix array file.
	protected static final String SUFFIX_EXTENSION = ".sarray";

	// the size of the bucket for the suffix array creation
	protected static final int BUCKET_SIZE = 5;

	// the size of an int primitive type in bytes
	protected static final int INT_BYTE_SIZE = Integer.SIZE / Byte.SIZE;



	/***** START OF TESTING AND DEBUGGING CODE, see here for examples of how to use the SuffixArray *****/
	/**
	 * Tester methods to test all substring can be retrieved. When testing make
	 * sure the sequence was created using a 1 to 1 mapping function of the
	 * alphabet to byte.
	 */
	private static void queryAllSubstrings(SuffixArray sa, SuffixArraySequence sequence, int iterations) {
		int tp = 0, fn = 0, tn = 0, fp = 0;

		Random r = new Random(); // random number generator
		for(int i = 0; i < iterations; i++) {
			int length = r.nextInt(50) + 5;  // random number from 5 to 40
			int position = r.nextInt((int)(sequence.getSize()-length));
			String query = sequence.getSubsequence(position, position+length);
			if(sequence.isEncodable(query)) {
				int pos = sa.search(sequence.toBytes(query));
				if(pos >= 0) {
					String match = sequence.getSubsequence(sa.getPosition(pos), sa.getPosition(pos)+length);
					if(match.equals(query)) {
						//System.out.println("We found correctly " + query + " at " + pos);
						//System.out.println(sequence.toString(sa.getPosition(pos), length));
						//System.exit(-1);
						tp++;  
					}
					else {
						fn++;
						System.out.println(query + '\t' + match);
					}
				}
				else {
					fn++;
					String match = sequence.getSubsequence(sa.getPosition(-pos-1), sa.getPosition(-pos-1)+length);
					System.out.println(query + "\t" + match);
				}
			}
			else {
				// nothing should be returned
				int pos = sa.search(sequence.toBytes(query));
				if(pos >= 0) {
					System.out.println("We found incorrectly " + query + " at " + pos);
					System.out.println(sequence.getSubsequence(sa.getPosition(pos), sa.getPosition(pos)+length));
					System.exit(-1);
					fp++;
				}
				else {
					tn++;
				}
			}
		}
		System.out.println();
		System.out.println("********** Test statistics **********");
		System.out.println("**** iterations: " + iterations);
		System.out.println("**** true positives: " + tp);
		System.out.println("**** false negative: " + fn);
		System.out.println("**** true negatives: " + tn);
		System.out.println("**** false positive: " + fp);
		System.out.println("**** sensitivity: " + (tp*100.0 / (tp + fn)));
		System.out.println("**** specificity: " + (tn*100.0 / (tn + fp)));
		System.out.println("*************************************");
		System.out.println();
	}


	/**
	 * Tester method.
	 */
	private static void debug() {
		String fastaFile;
		String userHome = System.getProperty("user.home");
		int iterations = 1000000;

		fastaFile = userHome+"/Data/Databases/tiny.fasta";
		//fastaFile = userHome+"/Data/Databases/small.fasta";
		//fastaFile = userHome+"/Data/Databases/single.fasta";
		//fastaFile = userHome+"/Data/Databases/uniprot_sprot.fasta";
		fastaFile = userHome+"/Data/Databases/yeast_nr050706.fasta";
		//fastaFile = "/home/sangtaekim/Research/Data/EColiDB/Ecol_protein_formatted.fasta";
		//fastaFile = "/home/sangtaekim/Research/Data/SProt/uniprot_sprot.fasta";
		//fastaFile = userHome+"/Desktop/test.fasta";
		//fastaFile = "/home/sangtaekim/Research/Data/HumanGenome/translated/HSRM.NCBI36.54.translation.0.fasta";

		long time = System.currentTimeMillis();
		SuffixArraySequence sequence = new SuffixArraySequence(fastaFile);
		System.out.println("-- Loading fasta file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");

		time = System.currentTimeMillis();
		SuffixArray sa = new SuffixArray(sequence);
		System.out.println("-- Loading SuffixArray file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");

		//MatchSet match = sa.findAll("PKVPFDPKFKEKLYDSYLDKAAKTK");

		//System.out.println("Translating _ " + sequence.toByte('_'));

		// print out the matches
		//for (int i=0; i< match.getSize(); i++) {
		//  System.out.println(sequence.toString(match.getStart(i), 10));
		//}

		time = System.currentTimeMillis();
		queryAllSubstrings(sa, sequence, iterations);

		System.out.println("-- Searching time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
	}



	/***** MEMBERS *****/
	// the indices of the sorted suffixes
	protected IntBuffer indices;

	// the sequence representing all the suffixes
	protected SuffixArraySequence sequence;

	// the class that generates suffixes from the given adapter
	protected SuffixFactory factory;

	// precomputed left-middle LCPs parameterized by the middle index
	protected ByteBuffer leftMiddleLcps;

	// precomputed middle-right LCPs parameterized by the middle index
	protected ByteBuffer middleRightLcps;

	// precomputed LCPs of neighboring suffixes
	protected ByteBuffer neighboringLcps;	// added by Sangtae

	// the number of suffixes in this array
	protected int size;



	/***** CLASS DEFINITION CODE *****/
	/**
	 * Print usage message.
	 */
	private static void printUsageAndExit() {
		System.out.println("usage: java SuffixArray [dbFile [queryFile]]");
		System.out.println("\tdbFile - the path to the database file with extension \".fasta\".");
		System.out.println("\tqueryFile - the path to the query file. One query per line. Use \"-\" for command line input.");
		System.out.println("\tArguments must be provided in order. Invocation with no arguments will run the tool through a series of test cases.");
		System.exit(-1);
	}


	/**
	 * Main method.
	 * @param args command line arguments.
	 */
	public static void main(String args[]) {

		if(args.length==0) {
			debug();
			return;
		}

		if(args.length<=2) {
			SuffixArraySequence sequence = new SuffixArraySequence(args[0]);
			SuffixArray sa = new SuffixArray(sequence);

			BufferedReader input = null;
			if(args.length==2) {
				try {
					input = new BufferedReader(new FileReader(args[1]));
				}
				catch(IOException e) {
					e.printStackTrace();
					System.exit(-1);
				}
			}
			else {
				input = new BufferedReader(new InputStreamReader(System.in));
			}

			sa.searchWithFile(input);
		}
		else {
			printUsageAndExit();
		}
	}


	/**
	 * Constructor that creates a suffixArray file from the given sequence. The
	 * name of the suffixArray will have the basePath of the sequence with the
	 * suffix array extension attached to it.
	 * @param sequence the sequence object to create the suffix array from.
	 * @param suffixFile the path to the precomputed suffix array file. If the
	 *        file does not exist, write it.
	 */
	public SuffixArray(SuffixArraySequence sequence, String suffixFile) {

		this.sequence = sequence;
		this.factory = new SuffixFactory(sequence);

		// create the file if it doesn't exist.
		if(!new File(suffixFile).exists()) {           
			createSuffixArrayFile(sequence, suffixFile);
		}

		// load the file
		int id = readSuffixArrayFile(suffixFile);

		// check that the files are consistent
		if(id != sequence.getId()) {
			System.err.println(suffixFile + " was not created from the sequence " + sequence.getBaseFilepath());
			System.err.println("Please recreate the suffix array file.");
			System.exit(-1);
		}

	}


	/**
	 * Constructor that attempts to read the suffix array from the provided file.
	 * @param sequence the sequence object.
	 */
	public SuffixArray(SuffixArraySequence sequence) {
		// infer the suffix array file from the sequence.
		this(sequence, sequence.getBaseFilepath() + SUFFIX_EXTENSION);
	}

	/** Constructor that reads the suffix array information from CompactSuffixArray
	 * @param
	 * @return
	 */
	public SuffixArray(CompactSuffixArray sa) {
		
	}
	
	
	public int getSize()
	{
		return size;
	}
	
	/**
	 * Helper function to initialize the leftMiddleLcps and middleRightLcps.
	 * @param nLcps the neigboring lcps.
	 * @param lLcps the left-middle lpcs.
	 * @param rLcps the middle=right lcps.
	 * @param start start index (inclusive).
	 * @param end end index (inclusive).
	 * @return the LCP between these two indices.
	 */
	private static byte initializeLcps(byte[] nLcps, byte[] lLcps, byte[] rLcps, int start, int end) {
		// base case
		if(end - start == 1) {
			// the assumption is that lcps[index] encodes the LCP(index-1, index)
			return nLcps[end]; 
		}

		// recursion
		int middleIndex = (start + end) / 2;
		byte lLcp = initializeLcps(nLcps, lLcps, rLcps, start, middleIndex);
		lLcps[middleIndex] = lLcp;
		byte rLcp = initializeLcps(nLcps, lLcps, rLcps, middleIndex, end);
		rLcps[middleIndex] = rLcp;

		// return the smallest one
		return lLcp < rLcp ? lLcp : rLcp;
	}



	/**
	 * Helper method that creates the suffixFile.
	 * @param sequence the Adapter object that represents the database (text).
	 * @param suffixFile the output file.
	 */
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
        }
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
			byte[] rLcps = new byte[(int)sequence.getSize()];
			byte[] lLcps = new byte[(int)sequence.getSize()];
			System.out.println("Computing the parameterized lcp arrays..");
			initializeLcps(neighboringLcps, lLcps, rLcps, 0, (int)(sequence.getSize()-1));
			out.write(lLcps);
			out.write(rLcps);
			out.write(neighboringLcps);	// Sangtae
			out.flush(); out.close();
		}
		catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return;
	}


	/**
	 * Helper method that initializes the suffixArray object from the file.
	 * Initializes indices, leftMiddleLcps, middleRightLcps and neighboringLcps.
	 * @param suffixFile the suffix array file.
	 * @return returns the id of this file for consistency check.
	 */
	protected int readSuffixArrayFile(String suffixFile) {
		try {
			// read the first integer which encodes for the size of the file
			DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(suffixFile)));
			this.size = in.readInt();
			// the second integer is the id
			int id = in.readInt();
			in.close();

			FileChannel fc = new FileInputStream(suffixFile).getChannel();

			// System.out.println("Reading the sorted indices.");
			long startPos = 2*INT_BYTE_SIZE;
			long sizeOfIndices = ((long)size)*INT_BYTE_SIZE;

			// read indices
			final int MAX_READ_SIZE = INT_BYTE_SIZE*(Integer.MAX_VALUE/4);
			IntBuffer[] dsts = new IntBuffer[(int)(sizeOfIndices/MAX_READ_SIZE)+1];
			for(int i=0; i<dsts.length; i++)
			{
				if(i<dsts.length-1)
				{
					dsts[i] = fc.map(FileChannel.MapMode.READ_ONLY, startPos, MAX_READ_SIZE).asIntBuffer();
					startPos += MAX_READ_SIZE;
				}
				else
				{
					dsts[i] = fc.map(FileChannel.MapMode.READ_ONLY, startPos, sizeOfIndices-(MAX_READ_SIZE)*(dsts.length-1)).asIntBuffer();
					startPos += sizeOfIndices-MAX_READ_SIZE*(dsts.length-1);
				}
			}

			if(dsts.length == 1)
				this.indices = dsts[0];
			else
			{
				// When sizeOfIndices > Integer.MAX_VALUE
				// It takes extra 5 seconds
				// totalCapacity must be smaller than Integer.MAX_VALUE
				long totalCapacity = 0;
				for(IntBuffer buf : dsts)
					totalCapacity += buf.capacity();
				assert(totalCapacity <= Integer.MAX_VALUE);
				//    	  System.out.println(totalCapacity);
				//   	  System.out.println(Runtime.getRuntime().totalMemory()+" " + Runtime.getRuntime().maxMemory()+" "+Runtime.getRuntime().freeMemory());
				this.indices = IntBuffer.allocate((int)totalCapacity);
				for(int i=0; i<dsts.length; i++)
				{
					for(int j=0; j<dsts[i].capacity(); j++)
						indices.put(dsts[i].get());
				}
				indices.rewind();
			}

			// System.out.println("Reading the leftMiddle lcps.");
			//      startPos += sizeOfIndices;
			int sizeOfLcps = size;
			this.leftMiddleLcps = fc.map(FileChannel.MapMode.READ_ONLY, startPos, sizeOfLcps).asReadOnlyBuffer();

			// System.out.println("Reading the middleRight lcps.");
			startPos += sizeOfLcps;
			this.middleRightLcps = fc.map(FileChannel.MapMode.READ_ONLY, startPos, sizeOfLcps).asReadOnlyBuffer();

			// added by Sangtae
			startPos += sizeOfLcps;
			this.neighboringLcps = fc.map(FileChannel.MapMode.READ_ONLY, startPos, sizeOfLcps).asReadOnlyBuffer();
			fc.close();

			return id;
		}
		catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

		return 0;
	}

	@Override
	public String toString() {
		String retVal = "Size of the suffix array: " + this.size + "\n";
		int rank = 0;
		while(indices.hasRemaining()) {
			int index = indices.get();
			int lcp = this.neighboringLcps.get(rank);
			retVal += rank + "\t" + index + "\t" + lcp + "\t" + sequence.toString(factory.makeSuffix(index).getSequence()) + "\n";
			rank++;
		}
		indices.rewind();        // reset marks after iteration
		neighboringLcps.rewind();
		return retVal;
	}


	/**
	 * This method translates the suffix array search index into a position of 
	 * the Adapter (sequence).
	 * @param index
	 */
	public int getPosition(int index) {
		if(index >= 0 && index < this.size)     return this.indices.get(index);
		return index;
	}


	/**
	 * Alternative to search the suffix array in which a MatchSet is return with
	 * all the starting positions in the sequence represented by this SuffixArray.
	 * @param pattern the ByteSequence to look for. The ByteSequence can be easily
	 *        translated from the Adapter sequence.
	 * @return a MatchSet object containing the match positions.
	 */
	public MatchSet findAll(ByteSequence pattern) {
		int matchIndex = search(pattern);
		MatchSet ms = new MatchSet();

		if(matchIndex >= 0) {
			for(int i = matchIndex; i < this.size; i++) {
				int start = getPosition(i);
				int numMatches = this.sequence.getLCP(pattern, start);
				if(numMatches == pattern.getSize())
					ms.add(start, start+pattern.getSize());
				else
					break;
			}
		}
		return ms;
	}

	/**
	 * Find all matches in the sequence represented by this SuffixArray and return their string representations.
	 * @author sangtaekim
	 * @param pattern the query string
	 * @return a list of matched strings.
	 */
	public ArrayList<String> getAllMatchedStrings(String pattern)
	{
		MatchSet matchSet = findAll(pattern);
		ArrayList<String> matches = new ArrayList<String>();
		for(int i=0; i<matchSet.getSize(); i++)
		{
			int start = matchSet.getStart(i);
			int end = matchSet.getEnd(i);
			matches.add(sequence.toChar(sequence.getByteAt(start-1))+"."+sequence.getSubsequence(start, end)+"."+sequence.toChar(sequence.getByteAt(end)));
		}
		return matches;
	}

	/**
	 * Find all matches in the sequence represented by this SuffixArray and return their string representations.
	 * @author sangtaekim
	 * @param pattern the query string
	 * @param lengthFlankingPep the length of flanking strings attached to the pattern
	 * @return a list of matched strings.
	 */
	public ArrayList<String> getAllMatchedStrings(String pattern, int lengthFlankingStr)
	{
		MatchSet matchSet = findAll(pattern);
		ArrayList<String> matches = new ArrayList<String>();
		for(int i=0; i<matchSet.getSize(); i++)
		{
			int start = matchSet.getStart(i);
			int end = matchSet.getEnd(i);
			String leftStr = sequence.getSubsequence(Math.max(0, start-lengthFlankingStr), start);
			String rightStr = sequence.getSubsequence(end+1, Math.min(end+1+lengthFlankingStr, sequence.getSize()));
			matches.add(leftStr+"."+sequence.getSubsequence(start, end)+"."+rightStr);
		}
		return matches;
	}

	/**
	 * Find all matches in the sequence represented by this SuffixArray and return annotations of all matched proteins.
	 * @param pattern the query string.
	 * @return a set of protein annotations.
	 */
	public ArrayList<String> getAllMatchingAnnotations(String pattern)
	{
		ArrayList<String> annotationSet = new ArrayList<String>();
		MatchSet matchSet = findAll(pattern);
		for(int i=0; i<matchSet.getSize(); i++)
			annotationSet.add(sequence.getAnnotation(matchSet.getStart(i)));

		return annotationSet;
	}

	/**
	 * Find the annotation of the corresponding index.
	 * @param pattern the query string.
	 * @return the annotation of the corresponding index.
	 */
	public String getAnnotation(int index)
	{
		return sequence.getAnnotation(index);
	}

	public ArrayList<String> getAllMatchingEntries(String pattern)
	{
		ArrayList<String> chunkSet = new ArrayList<String>();
		MatchSet matchSet = findAll(pattern);
		for(int i=0; i<matchSet.getSize(); i++)
			chunkSet.add(sequence.getMatchingEntry(matchSet.getStart(i)));

		return chunkSet;
	}

	/**
	 * 
	 * @param pattern
	 * @return
	 */
	public MatchSet findAll(String pattern) {
		return findAll(sequence.toBytes(pattern));
	}


	/**
	 * Alternative method of searching that takes input as a string.
	 * @param pattern the pattern in String form.
	 * @return the index returned is the relative position in this suffix array. To
	 * get the index in the Adapter sequence, call getPosition.
	 */
	public int search(String pattern) {
		return search(sequence.toBytes(pattern));
	}

	/**
	 * <p>The generalized search method for this suffixArray. This search routine
	 * does a binary search on the suffixArray and returns the starting index
	 * of the pattern. A positive number indicates a successful match, while a
	 * negative return value means no match.</p>
	 * <p>It is very easy to decode the match indices. The return value is
	 * guaranteed to be the left-most (smallest) match in the suffix array. 
	 * Therefore, to retrieve all the matches, one only needs to walk to the right
	 * until the sorted suffixes do not match the query.</p>
	 * <p>For negative values, it represents the insertion point of the pattern
	 * into the suffix array shifted by 1. For example, if the return value is m,
	 * then pattern should be inserted at -m-1 and all elements at -m-1, including
	 * -m-1 shifted to the right by 1 position. In other words element at -m-1 is
	 * the first element that is lexographically greater that pattern.</p>
	 * <p>This implementation takes O(P+logN) per execution, where P is the length
	 * of the pattern and N is the size of the suffix array (Manber & Myers method).</p> 
	 * @param pattern the query to search for in the suffix array.
	 * @return the index returned is the relative position in this suffix array. To
	 * get the index in the Adapter sequence, call getPosition.
	 */
	public int search(ByteSequence pattern) {

		// check that the pattern is within the left boundary
		int leftResult = pattern.compareTo(factory.makeSuffix(indices.get(0)));
		if(Math.abs(leftResult)-1 == pattern.getSize())
			return 0;              // exact leftmost match of the first element
		if(leftResult < 0)       
			return -1;             // insertion point is at position 0

		// check that the pattern is within the right boundary
		int rightResult = factory.makeSuffix(indices.get(this.size-1)).compareTo(pattern); 
		if(rightResult < 0)      
			return -this.size;     // insertion point is at the end of the array

		// initialize the longest common prefixes values
		int queryLeftLcp = pattern.getLCP(factory.makeSuffix(indices.get(0)));
		int queryRightLcp = pattern.getLCP(factory.makeSuffix(indices.get(this.size-1)));

		// debug code
		// System.out.println(queryLeftLcp + "\t" + queryRightLcp);

		// indices for the binary search
		int leftIndex = 0;
		int rightIndex = (int)sequence.getSize()-1;

		// loop invariant: element at leftIndex < pattern <= element at rightIndex
		while(rightIndex-leftIndex > 1) {

			int middleIndex = (leftIndex+rightIndex)/2;
			if(queryLeftLcp >= queryRightLcp) {
				byte leftMiddleLcp = this.leftMiddleLcps.get(middleIndex);
				if(leftMiddleLcp > queryLeftLcp) {       // and queryMiddle == queryLeft
					leftIndex = middleIndex;
					// queryLeft = queryMiddle, already true
				}
				else if(queryLeftLcp > leftMiddleLcp) {  // and queryMiddle == leftMiddle
					// we can conclude that query < middle because queryMiddle < queryLeft
					queryRightLcp = leftMiddleLcp;
					rightIndex = middleIndex;
				}
				else {     // queryLeft == leftMiddle == queryMiddle 
					int middleResult = Math.min(pattern.compareTo(factory.makeSuffix(indices.get(middleIndex)), queryLeftLcp), Byte.MAX_VALUE);
					if(middleResult <= 0) {      // pattern <= middle
						queryRightLcp = middleResult == 0 ? pattern.getSize() : -middleResult-1;
						rightIndex = middleIndex;
					}
					else {                       // middle < pattern
						queryLeftLcp = middleResult-1;
						leftIndex = middleIndex;
					}
				}
			}
			else {       // queryRight > queryLeft
				int middleRightLcp = this.middleRightLcps.get(middleIndex);
				if(middleRightLcp > queryRightLcp) {           // and queryMiddle == queryRight
					rightIndex = middleIndex;
					// queryRight = queryMiddle, already true
				}
				else if(queryRightLcp > middleRightLcp) {      // and queryMiddle == middleRight
					queryLeftLcp = middleRightLcp;
					leftIndex = middleIndex;
				}
				else {     // middleRight == queryRight == queryMiddle
					int middleResult = Math.min(pattern.compareTo(factory.makeSuffix(indices.get(middleIndex)), queryRightLcp), Byte.MAX_VALUE);
					if(middleResult <= 0) {      // pattern <= middle
						queryRightLcp = middleResult == 0 ? pattern.getSize() : -middleResult-1;
						rightIndex = middleIndex;
					}
					else {                       // middle < pattern
						queryLeftLcp = middleResult-1;
						leftIndex = middleIndex;
					}
				}
			}
		}

		// evaluate the base cases, found!
		if(queryRightLcp == pattern.getSize())   return rightIndex;

		// not found
		return -rightIndex-1;
	}


	// To-do (sangtae): search suffix array using partial matches
	// public MatchSet search(MatchSet partialMatch, byte b)

	/**
	 * Treat the parameter as the source of input. One line per query.
	 * @param in the queries. One input per line. IMPLEMENT!!!
	 */
	public void searchWithFile(BufferedReader in) {
		return;
	}

	public void printAllPeptides(AminoAcidSet aaSet, int minLength, int maxLength)
	{
		//		  ArrayList<Pair<Float,Integer>> pepList = new ArrayList<Pair<Float,Integer>>();

		double[] aaMass = new double[128];
		for(int i=0; i<aaMass.length; i++)
			aaMass[i] = -1;
		for(AminoAcid aa : aaSet)
			aaMass[aa.getResidue()] = aa.getAccurateMass();
		double[] prm = new double[maxLength];
		int rank = 0;
		int i = Integer.MAX_VALUE;
		while(indices.hasRemaining()) {
			int index = indices.get();
			int lcp = this.neighboringLcps.get(rank);
			rank++;
			//			  System.out.println(sequence.getSubsequence(index, index+10)+":"+index+":"+lcp);
			if(lcp > i)
				continue;
			for(i=lcp; i<maxLength; i++)
			{
				char residue = sequence.getCharAt(index+i);
				double m = aaMass[residue];
				if(m <= 0)
				{
					break;
				}
				if(i != 0)
					prm[i] = prm[i-1] + m;
				else
					prm[i] = m;
				if(i+1 >= minLength && i+1 <= maxLength)
					//					  ;
				//					  pepList.add(new Pair<Float,Integer>((float)prm[i], index));
					System.out.println(index+"\t"+(float)prm[i]+"\t"+sequence.getSubsequence(index, index+i+1));
			}
		}

		//		  Collections.sort(pepList, new Pair.PairComparator<Float,Integer>());
		//		  System.out.println("Sorted");
		indices.rewind();
		neighboringLcps.rewind();
	}



	public int getNumCandidatePeptides(AminoAcidSet aaSet, float peptideMass, Tolerance tolerance)
	{
		double[] aaMass = new double[128];
		for(int i=0; i<aaMass.length; i++)
			aaMass[i] = -1;
		for(AminoAcid aa : aaSet)
			aaMass[aa.getResidue()] = aa.getAccurateMass();
		int maxLength = 50;
		float tolDa = tolerance.getToleranceAsDa(peptideMass);
		double[] prm = new double[maxLength];
		int numCandidatePeptides = 0;
		int rank = 0;
		int matchLength = Integer.MAX_VALUE;
		while(indices.hasRemaining()) {
			int index = indices.get();
			int lcp = this.neighboringLcps.get(rank);
			//			  System.out.println(sequence.getSubsequence(index, index+10)+":"+index+":"+lcp);
			rank++;
			if(lcp >= matchLength)
			{
				numCandidatePeptides++;
				continue;
			}
			for(int i=lcp; i<maxLength; i++)
			{
				char residue = sequence.getCharAt(index+i);
				double m = aaMass[residue];
				if(m <= 0)
				{
					matchLength = Integer.MAX_VALUE;
					break;
				}
				//				  if(sequence.getSubsequence(index, index+10).contains("_"))
				//					  System.out.println("Debug");
				if(i != 0)
					prm[i] = prm[i-1] + m;
				else
					prm[i] = m;
				if(prm[i] <= peptideMass - tolDa)
					continue;
				else if(prm[i] < peptideMass + tolDa)
				{
					matchLength = i;
					numCandidatePeptides++;
					break;
				}
				else
				{
					matchLength = Integer.MAX_VALUE;
					break;
				}
			}
		}

		indices.rewind();
		neighboringLcps.rewind();
		return numCandidatePeptides;
	}		
	
	/***** METHODS NOT PORTED *****
	 *  
  public boolean canThrowOut(String seq) {
    return search(seq) < 0;
  }


  public String searchForString(String pattern)
  {
    int index = search(pattern);
    if(index < 0)
      return null;
    else
    {
      int startPos = pos[index];
      return getMatchedString(startPos, pattern.length());
    }
  }

  public String[] searchForAllStringWithAnnotation(String pattern)
  {
    int index = search(pattern);
    if(index < 0)
      return null;
    else
    {
      TreeMap<Integer,String> matches = new TreeMap<Integer,String>();
      int startPos = pos[index];
      int minPos = startPos;
      int maxPos = startPos;
      String match = getMatchedString(startPos, pattern.length());
      String peptide = match.substring(match.indexOf('.')+1,match.lastIndexOf('.'));
      matches.put(startPos, match+"\t"+getAnnotation(startPos));

      for(int i=index-1; i>=0; i--)
      {
        startPos = pos[i];
        String matchedString = getMatchedString(startPos, pattern.length());
        String matchedPeptide = matchedString.substring(matchedString.indexOf('.')+1, matchedString.lastIndexOf('.'));
        boolean isMatch = true;
        for(int l=pattern.length()-1; l>=0; l--)
        {
          if(HashIndex.getAAIndex(matchedPeptide.charAt(l)) != HashIndex.getAAIndex(pattern.charAt(l)))
          {
            isMatch = false;
            break;
          }
        }
        if(isMatch)
        {
          matches.put(startPos, matchedString+"\t"+getAnnotation(startPos));
        }
        else
          break;
      }
      for(int i=index+1; i+pattern.length()<alphabet.length; i++)
      {
        startPos = pos[i];
        String matchedString = getMatchedString(startPos, pattern.length());
        String matchedPeptide = matchedString.substring(matchedString.indexOf('.')+1, matchedString.lastIndexOf('.'));
        boolean isMatch = true;
        for(int l=pattern.length()-1; l>=0; l--)
        {
          if(HashIndex.getAAIndex(matchedPeptide.charAt(l)) != HashIndex.getAAIndex(pattern.charAt(l)))
          {
            isMatch = false;
            break;
          }
        }
        if(isMatch), tokens[1]
          matches.put(startPos, matchedString+"\t"+getAnnotation(startPos));
        else
          break;
      }
      ArrayList<String> results = new ArrayList<String>();
      Iterator<Map.Entry<Integer, String>> itr = matches.entrySet().iterator();
      while(itr.hasNext())
        results.add(itr.next().getValue());
      return results.toArray(new String[0]);
    }
  }

  public String searchForStringWithAnnotation(String pattern)
  {
    int index = search(pattern);
    if(index < 0)
      return null;
    else
    {
      int startPos = pos[index];

      return getMatchedString(startPos, pattern.length())+"\t"+getAnnotation(startPos);
    }
  }

  public String getAnnotationByIndex(int index, int length)
  {
    int startPos = pos[index];
    return getMatchedString(startPos, length);
  }

  public String getMatchedString(int startPos, int length)
  {
    StringBuffer str = new StringBuffer();
    if(startPos > 0)
    {
      if(origSeq[startPos-1] >= 0)
        str.append(HashIndex.getAAFromIndex20(origSeq[startPos-1]));
    }
    str.append(".");
    for(int i=startPos; i<startPos+length; i++)
      str.append((HashIndex.getAAFromIndex20(origSeq[i])));
    str.append(".");
    if(startPos+length < origSeq.length)
    {
      if(origSeq[startPos+length] >= 0)
        str.append(HashIndex.getAAFromIndex20(origSeq[startPos+length]));
    }
    return str.toString();
  }

  public String getAnnotation(int startPos)
  {
    String annotation = annotations.floorEntry(startPos).getValue();
    return annotation;
  }

	public static void indexFastaFile(String fileName)
	{
		String name = fileName.substring(0, fileName.lastIndexOf('.'));
		serializeFasta(fileName, name+".serial", name+".annotation");
		generateSuffixArray(name+".serial", name+".sarr", name+".lcp");
	}

	public static void indexTextFile(String fileName)
	{
		String name = fileName.substring(0, fileName.lastIndexOf('.'));
		serializeSequences(fileName, name+".serial");
		generateSuffixArray(name+".serial", name+".sarr", name+".lcp");
	}
	 */

}
