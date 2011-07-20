package msdbsearch;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.IntBuffer;
import java.nio.channels.FileChannel;

import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;

public class SuffixArrayForMSGFDB extends SuffixArray {

	private int[] numDisinctPeptides;
	
	public SuffixArrayForMSGFDB(SuffixArraySequence sequence, int minPeptideLength, int maxPeptideLength) {
		super(sequence);
		
		// compute the number of distinct peptides
		numDisinctPeptides = new int[maxPeptideLength+2];
		for(int length=minPeptideLength; length<=maxPeptideLength+1; length++)
			numDisinctPeptides[length] = getNumDistinctSeq(length);
	}

	public IntBuffer getIndices()	{	return indices; }
	public ByteBuffer getNeighboringLcps()	{	return neighboringLcps; }
	public SuffixArraySequence getSequence()	{ return sequence; }
	public int getNumDistinctPeptides(int length)	{ return numDisinctPeptides[length]; }
	
	@Override
	protected int readSuffixArrayFile(String suffixFile) {
		try {
			// read the first integer which encodes for the size of the file
			DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(suffixFile)));
			size = in.readInt();
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

			int sizeOfLcps = size;
			// leftMiddleLcps are not read.
//			this.leftMiddleLcps = fc.map(FileChannel.MapMode.READ_ONLY, startPos, sizeOfLcps).asReadOnlyBuffer();

			startPos += sizeOfLcps;
			// middleRightLcps are not read.
//			this.middleRightLcps = fc.map(FileChannel.MapMode.READ_ONLY, startPos, sizeOfLcps).asReadOnlyBuffer();

			// added by Sangtae
			startPos += sizeOfLcps;
			neighboringLcps = fc.map(FileChannel.MapMode.READ_ONLY, startPos, sizeOfLcps).asReadOnlyBuffer();
//			neighboringLcps = 
			fc.close();

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
