package edu.ucsd.msjava.msgf;

public class IntHistogram extends Histogram<Integer> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	// assuming the hisgram is centered around zero
	public float[] getSmoothedHist(int keySize)
	{
		float[] smoothedHist = new float[keySize*2+1];
		// smoothing
		for(int i=-keySize; i<=keySize; i++)
		{
			int windowSize;
			if(Math.abs(i)<=3)
				windowSize = 0;
			else
				windowSize = (int)(Math.log(Math.abs(i))/Math.log(2))-1;
			
			int numUsedEntries = 0;
			int sum = 0;
			
			for(int j=i-windowSize; j<=i+windowSize; j++)
			{
				if(j<=keySize && j>=-keySize)
				{
					numUsedEntries++;
					sum += this.get(j);
				}
			}

			while(sum == 0)
			{
				windowSize++;
				if(windowSize > keySize)
				{
					sum = 1;
					numUsedEntries = 2*keySize+1;
					break;
				}
				else
				{
					if(i-windowSize >= -keySize)
					{
						sum += this.get(i-windowSize);
						numUsedEntries++;
					}
					if(i+windowSize <= keySize)
					{
						sum += this.get(i+windowSize);
						numUsedEntries++;
					}
				}
			}
			smoothedHist[i+keySize] = sum/(float)numUsedEntries;
		}
		return smoothedHist;
	}

}
