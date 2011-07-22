package msdbsearch;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.IntBuffer;

import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;

public class SuffixArrayForMSGFDB extends SuffixArray {

	private int[] numDisinctPeptides;
	
	public SuffixArrayForMSGFDB(SuffixArraySequence sequence)
	{
		super(sequence);
	}
	
	public SuffixArrayForMSGFDB(SuffixArraySequence sequence, int minPeptideLength, int maxPeptideLength) 
	{
		super(sequence);
		
		// compute the number of distinct peptides
		numDisinctPeptides = new int[maxPeptideLength+2];
		for(int length=minPeptideLength; length<=maxPeptideLength+1; length++)
			numDisinctPeptides[length] = getNumDistinctSeq(length);
	}

	public IntBuffer getIndices()	{	return indices; }
	public ByteBuffer getNeighboringLcps()	{	return neighboringLcps; }
	public SuffixArraySequence getSequence()	{ return sequence; }
	public int getNumDistinctPeptides(int length)	
	{
		if(numDisinctPeptides != null)
			return numDisinctPeptides[length];
		else
			return this.getNumDistinctSeq(length);
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
			indices = IntBuffer.wrap(indexArr);

			int sizeOfLcps = size;
			// skip leftMiddleLcps and middleRightLcps
			in.skip(2*sizeOfLcps);
			// read neighboringLcps
			byte[] neighboringLcpArr = new byte[sizeOfLcps];
			in.read(neighboringLcpArr);
			neighboringLcps = ByteBuffer.wrap(neighboringLcpArr);
			in.close();

			return id;
		}
		catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

		return 0;
	}		
	
	private int getNumDistinctSeq(int length)
	{
		int numDistinctSeq = 0;
		while(neighboringLcps.hasRemaining())
		{
			int lcp = neighboringLcps.get();
			if(lcp < length)
			{
				numDistinctSeq++;
			}
		}
		neighboringLcps.rewind();
		indices.rewind();
		return numDistinctSeq++;
	}
}
