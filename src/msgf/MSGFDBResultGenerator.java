package msgf;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

public class MSGFDBResultGenerator extends ArrayList<MSGFDBResultGenerator.DBMatch> {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String header;
	public MSGFDBResultGenerator(String header)
	{
		this.header = header;
	}
	public void computeEFDR()
	{
		Collections.sort(this);
		float eTD = 0;	// expected target discovery
		float cumulativePValue = 0;
		for(int i=0; i<this.size(); i++)
		{
			float specProb = get(i).getSpecProb();
			cumulativePValue += get(i).getPValue();
			float eDD = cumulativePValue;	// expected decoy discovery
			for(int j=i+1; j<this.size(); j++)
				eDD += get(j).getEDD(specProb);
			eTD += 1-get(i).getPValue();
			get(i).setEFDR(Math.min(eDD/eTD, 1));
		}
	}
	
	public void writeResults(PrintStream out)
	{
		out.println(header+"\tEFDR");
		for(MSGFDBResultGenerator.DBMatch m : this)
			out.println(m.getResultStr()+"\t"+m.getEFDR());
	}
	
	public static class DBMatch implements Comparable<DBMatch>
	{
		private float specProb;
		private float pValue;
		private int numPeptides;
		private String resultStr;
		private float[] cumScoreDist;
		private float eFDR;
		
		public DBMatch(float specProb, int numPeptides, String resultStr, ScoreDist scoreDist) {
			this.specProb = specProb;
			this.pValue = getPValue(specProb, numPeptides);
			this.numPeptides = numPeptides;
			this.resultStr = resultStr;
			
			if(scoreDist.isProbSet())
			{
				this.cumScoreDist = new float[scoreDist.getMaxScore()-scoreDist.getMinScore()];
				cumScoreDist[0] = scoreDist.getProbability(scoreDist.getMaxScore()-1);
				int index = 1;
				for(int t=scoreDist.getMaxScore()-2; t>=scoreDist.getMinScore(); t--)
				{
					cumScoreDist[index] = cumScoreDist[index-1] + scoreDist.getProbability(t);
					index++;
				}
			}
		}
		
		public static float getPValue(float specProb, int numPeptides)
		{
			float pValue;
			double probCorr = 1.-(double)specProb;
			if(probCorr < 1.)
				pValue = (float)(1.- Math.pow(probCorr, numPeptides));
			else
				pValue = specProb*numPeptides;
			return pValue;
		}
		
		public void setEFDR(float eFDR)	{ this.eFDR = eFDR; }
		
		public float getEFDR() {
			return eFDR;
		}
		
		/**
		 * Gets expected decoy discovery for a given specProbThreshold
		 */
		public float getEDD(float specProbThreshold)	{
			float probEqualOrBetterTargetPep;
			if(specProbThreshold >= specProb)
				probEqualOrBetterTargetPep = specProb;
			else
				probEqualOrBetterTargetPep = getSpectralProbability(specProbThreshold);
			
			float pValue = getPValue(probEqualOrBetterTargetPep, numPeptides);
			return pValue;
		}
		
		// returns cumulative probability <= specProbThreshold
		public float getSpectralProbability(float specProbThreshold)
		{
			int index = Arrays.binarySearch(cumScoreDist, specProbThreshold);
			if(index >= 0)
				return cumScoreDist[index];
			else
			{
				index = -index-1;
				if(index > 0)
					return cumScoreDist[index-1];
				else
					return 0;
			}
		}
		
		public float getSpecProb() {
			return specProb;
		}
		public float getPValue() {
			return pValue;
		}
		public String getResultStr() {
			return resultStr;
		}
		@Override
		public int compareTo(DBMatch arg0) {
			if(this.specProb < arg0.specProb)
				return -1;
			else if(this.specProb > arg0.specProb)
				return 1;
			else
				return 0;
		}
	}
}
