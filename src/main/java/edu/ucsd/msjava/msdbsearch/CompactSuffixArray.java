package edu.ucsd.msjava.msdbsearch;

import java.io.*;
import java.util.*;

import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.sequences.Constants;
import edu.ucsd.msjava.suffixarray.SuffixFactory;

/**
 * SuffixArray class for fast exact matching.
 * @author Sangtae Kim
 *
 */
public class CompactSuffixArray {

	public static final int COMPACT_SUFFIX_ARRAY_FILE_FORMAT_ID = 8258;
	
	/***** CONSTANTS *****/
	// The default extension of a suffix array file.
	protected static final String EXTENSION_INDICES = ".csarr";
	protected static final String EXTENSION_NLCPS = ".cnlcp";

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
	private int maxPeptideLength;
	private int[] numDisinctPeptides;
	
	
	/**
	 * Constructor that attempts to read the suffix array from the provided file.
	 * @param sequence the sequence object.
	 */
	public CompactSuffixArray(CompactFastaSequence sequence) {
		// infer the suffix array file from the sequence.
		this.sequence = sequence;
		this.size = (int)sequence.getSize();
		this.factory = new SuffixFactory(sequence);
		indexFile = new File(sequence.getBaseFilepath() + EXTENSION_INDICES);
		nlcpFile = new File(sequence.getBaseFilepath() + EXTENSION_NLCPS);

		// create the file if it doesn't exist.
		if(!indexFile.exists() || !nlcpFile.exists() || !isCompactSuffixArrayValid()) {           
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
	}
	
	/**
	 * Constructor that attempts to read the suffix array from the provided file.
	 * @param sequence the sequence object.
	 */
	public CompactSuffixArray(CompactFastaSequence sequence, int maxPeptideLength) {
		this(sequence);
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
	
	private boolean isCompactSuffixArrayValid()
	{
		File[] files = {indexFile, nlcpFile};
		
		for(File f : files)
		{
			try {
				RandomAccessFile raf = new RandomAccessFile(f, "r");
				raf.seek(raf.length()-Integer.SIZE/8);
				int id = raf.readInt();
				raf.close();
				
				if(id != COMPACT_SUFFIX_ARRAY_FILE_FORMAT_ID)
					return false;
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		return true;
	}
	
	//TODO: this method has a bug
	private void computeNumDistinctPeptides()
	{
		boolean[] isValidResidue = new boolean[128];
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSet();
		for(AminoAcid aa : aaSet)
			isValidResidue[aa.getResidue()] = true;

		numDisinctPeptides = new int[maxPeptideLength+2];
		try {
			DataInputStream indices = new DataInputStream(new BufferedInputStream(new FileInputStream(getIndexFile())));
			indices.skip(CompactSuffixArray.INT_BYTE_SIZE*2);	// skip size and id
			
			DataInputStream neighboringLcps = new DataInputStream(new BufferedInputStream(new FileInputStream(nlcpFile)));
			int size = neighboringLcps.readInt();
			neighboringLcps.readInt();	// skip id
			
			for(int i=0; i<size; i++)
			{
				int index = indices.readInt();
				byte lcp = neighboringLcps.readByte();
				int idx = sequence.getCharAt(index); 
				if(isValidResidue[idx] == false)
					continue;

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
//			private static final int INCREMENT_SIZE = 10;
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
//					this.items = Arrays.copyOf(this.items, this.size+INCREMENT_SIZE);
					this.items = Arrays.copyOf(this.items, this.size*2);
//					if(this.size > 100)
//						System.out.println(item+":"+this.size);
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
		System.out.println("AlphabetSize: " + sequence.getAlphabetSize());
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
			if(j % 10000001 == 0)  System.out.printf("Suffix creation: %.2f%% complete.\n", j*100.0/sequence.getSize());

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

		System.gc();
		
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
				if(i % 1000000 == 0)   System.out.printf("Sorting %.2f%% complete.\n", i*100.0/bucketSuffixes.length);

				if(bucketSuffixes[i] != null) {

					SuffixFactory.Suffix[] sortedSuffixes = bucketSuffixes[i].getSortedSuffixes();

					SuffixFactory.Suffix first = sortedSuffixes[0];
					byte lcp = 0;
					if(prevBucketSuffix != null) {
						lcp = first.getLCP(prevBucketSuffix);
					}
					// write information to file
					indexOut.writeInt(first.getIndex());
					nLcpOut.writeByte(lcp);
					SuffixFactory.Suffix prevSuffix = first;

					for(int j = 1; j < sortedSuffixes.length; j++) {
						SuffixFactory.Suffix thisSuffix = sortedSuffixes[j];
						//store the information
						indexOut.writeInt(thisSuffix.getIndex());
						lcp = thisSuffix.getLCP(prevSuffix, BUCKET_SIZE);
						nLcpOut.writeByte(lcp);
						prevSuffix = thisSuffix;
					}
					prevBucketSuffix = sortedSuffixes[0];
					bucketSuffixes[i] = null;	// deallocate the memory
				}
			}

			bucketSuffixes = null;
			
			indexOut.writeInt(CompactSuffixArray.COMPACT_SUFFIX_ARRAY_FILE_FORMAT_ID);
			indexOut.flush();
			indexOut.close();
			
			nLcpOut.writeInt(CompactSuffixArray.COMPACT_SUFFIX_ARRAY_FILE_FORMAT_ID);
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
	
	public void measureNominalMassError(AminoAcidSet aaSet) throws Exception
	{
		//		  ArrayList<Pair<Float,Integer>> pepList = new ArrayList<Pair<Float,Integer>>();
		double[] aaMass = new double[128];
		int[] nominalAAMass = new int[128];
		for(int i=0; i<aaMass.length; i++)
		{
			aaMass[i] = -1;
			nominalAAMass[i] = -1;
		}
		
		for(AminoAcid aa : aaSet)
		{
			aaMass[aa.getResidue()] = aa.getAccurateMass();
			nominalAAMass[aa.getResidue()] = aa.getNominalMass();
		}
		double[] prm = new double[maxPeptideLength];
		int[] nominalPRM = new int[maxPeptideLength];
		int i = Integer.MAX_VALUE-1000;
		int[] numPeptides = new int[maxPeptideLength];
		int[][] numPepWithError = new int[maxPeptideLength][11];
		
		DataInputStream indices = new DataInputStream(new BufferedInputStream(new FileInputStream(getIndexFile())));
		indices.skip(CompactSuffixArray.INT_BYTE_SIZE*2);	// skip size and id

		DataInputStream nlcps = new DataInputStream(new BufferedInputStream(new FileInputStream(getNeighboringLcpFile())));
		nlcps.skip(CompactSuffixArray.INT_BYTE_SIZE*2);

		int size = this.getSize();
		int index = -1;
		for(int bufferIndex=0; bufferIndex<size; bufferIndex++) {
			index = indices.readInt();
			int lcp = nlcps.readByte();
			
			int idx = sequence.getCharAt(index); 
			if(aaMass[idx] <= 0)
				continue;
			
			if(lcp > i)
				continue;
			for(i=lcp; i<maxPeptideLength; i++)
			{
				char residue = sequence.getCharAt(index+i);
				double m = aaMass[residue];
				if(m <= 0)
				{
					break;
				}
				if(i != 0)
				{
					prm[i] = prm[i-1] + m;
					nominalPRM[i] = nominalPRM[i-1] + nominalAAMass[residue];
				}
				else
				{
					prm[i] = m;
					nominalPRM[i] = nominalAAMass[residue];
				}
				if(i+1 <= maxPeptideLength)
				{
					numPeptides[i]++;
					int error = (int)Math.round(prm[i]*0.9995)-nominalPRM[i];
					error += 5;
					numPepWithError[i][error]++;
//					System.out.println(index+"\t"+(float)prm[i]+"\t"+sequence.getSubsequence(index, index+i+1));
				}
			}
		}

		long total = 0;
		long totalErr = 0;
		System.out.println("Length\tNumDistinctPeptides\tNumPeptides\tNumPeptidesWithErrors");
		for(i=0; i<maxPeptideLength; i++)
		{
			System.out.print((i+1)+"\t"+this.numDisinctPeptides[i+1]+"\t"+numPeptides[i]);
			total += numPeptides[i];
			for(int j=0; j<11; j++)
			{
				if(numPepWithError[i][j]>0)
				{
					System.out.print("\t"+(j-5)+":"+numPepWithError[i][j]);
					if(j != 5)
						totalErr += numPepWithError[i][j];
				}
			}
			System.out.println("\t"+total+"\t"+totalErr+"\t"+(totalErr/(double)total));
		}
		System.out.println("Total #Peptides\t" + total);
		System.out.println("Total #Peptides with nominalMass errors\t" + totalErr + "\t" + totalErr/(double)total);
		
		indices.close();
		nlcps.close();
	}		
}
