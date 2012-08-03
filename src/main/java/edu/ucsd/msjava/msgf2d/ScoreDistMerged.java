package edu.ucsd.msjava.msgf2d;

public class ScoreDistMerged {
	private ScoreDist2D dist2D;
	
	// #better scoring peptides at score (t1, t2) where t1 >= i && t2 >= j
	private float[][] numPeptidesWithBetterScores;
	// total probability of peptides at score (t1, t2) where t1 >= i && t2 >= j
	private float[][] specProbWithBetterScores;
	
	// #better scoring peptides where score1 >= t1 && score2 == t2
	private float[][] numPeptidesBetterScore1;
	private float[][] probBetterScore1;
	
	// spectral probability
	private float[] specProb1;
	private float[] specProb2;
	
	// #peptides with equal or better scores
	@SuppressWarnings("unused")
	private float[] numBetter1;
	@SuppressWarnings("unused")
	private float[] numBetter2;
	
	// #better scoring peptides where score1 == t1 && score2 >= t2
	private float[][] numPeptidesBetterScore2;
	private float[][] probBetterScore2;
	
	private int minScore1, minScore2;
	private int maxScore1, maxScore2;
	
	public ScoreDistMerged(ScoreDist2D dist2D)
	{
		this.dist2D = dist2D;
		minScore1 = dist2D.scoreBound1.getMinScore();
		minScore2 = dist2D.scoreBound2.getMinScore();
		maxScore1 = dist2D.scoreBound1.getMaxScore()-1;
		maxScore2 = dist2D.scoreBound2.getMaxScore()-1;
		initBothBetter();
	}

	public int getMaxScore1()	{ return dist2D.getMaxScore1(); }
	public int getMaxScore2()	{ return dist2D.getMaxScore2(); }
	
	public int getMinScore1()	{ return dist2D.getMinScore1(); }
	public int getMinScore2()	{ return dist2D.getMinScore2(); }

	public float getProbabilityAt(int score1, int score2)	{ return dist2D.getProbability(score1, score2); }
	public float getNumRecsAt(int score1, int score2)	{ return dist2D.getNumRecs(score1, score2); }
	
	// t1 >= score1, t2 >= score2
	public float getNumBetterBoth(int score1, int score2)
	{
		return numPeptidesWithBetterScores[score1-minScore1][score2-minScore2];
	}
	
	public float getProbBetterBoth(int score1, int score2)
	{
		return specProbWithBetterScores[score1-minScore1][score2-minScore2];
	}
	
	// t1 >= score1
	public float getNumEqualOrBetterPeptides1(int score1)
	{
		int index1 = score1 - minScore1;
		float numRecs = 0;
		for(int score2=dist2D.getMinScore2(); score2<dist2D.getMaxScore2(); score2++)
			numRecs += numPeptidesBetterScore1[index1][score2-minScore2];
		return numRecs;
	}

	public float getSpectralProbability1(int score1)
	{
		int index1 = score1 - minScore1;
		
		if(score1 < minScore1)
			return 1;
		else if(score1 > maxScore1)
			return 0;
		else
			return specProb1[index1];
	}

	// t2 >= score2
	public float getNumEqualOrBetterPeptides2(int score2)
	{
		int index2 = score2 - minScore2;
		float numRecs = 0;
		for(int score1=dist2D.getMinScore1(); score1<dist2D.getMaxScore1(); score1++)
			numRecs += numPeptidesBetterScore2[score1-minScore1][index2];
		return numRecs;
	}

	public float getSpectralProbability2(int score2)
	{
		int index2 = score2 - minScore2;
		if(score2 < minScore2)
			return 1;
		else if(score2 > maxScore2)
			return 0;
		else
			return specProb2[index2];
	}

	// use specProb(t1)*specProb(t2) as ranking method
	public float getNumEqualOrBetterPeptides(int score1, int score2)
	{
		float numBetter = 0;

		int index1 = score1-minScore1;
		int index2 = score2-minScore2;
		
		numBetter += numPeptidesWithBetterScores[index1][index2];
		float curCombinedSpecProb = specProb1[index1]*specProb2[index2];
		
		// t1 < score1 && t2 > score2
		int t2 = maxScore2;
		for(int t1=minScore1; t1<score1; t1++)
		{
			for(; t2>=score2; t2--)
			{
				float combinedSpecProb = specProb1[t1-minScore1]*specProb2[t2-minScore2];
				if(combinedSpecProb > curCombinedSpecProb)
					break;
			}
			if(t2+1 <= maxScore2)
				numBetter += numPeptidesBetterScore2[t1-minScore1][t2+1-minScore2];
		}
		
		// t1 > score1 && t2 < score2
		int t1 = maxScore1;
		for(t2=minScore2; t2<score2; t2++)
		{
			for(; t1>=score1; t1--)
			{
				float combinedSpecProb = specProb1[t1-minScore1]*specProb2[t2-minScore2];
				if(combinedSpecProb > curCombinedSpecProb)
					break;
			}
			if(t1+1 <= maxScore1)
				numBetter  += numPeptidesBetterScore1[t1+1-minScore1][t2-minScore2];
		}		
		return numBetter;
	}	
	
	// use specProb(t1)*specProb(t2) as ranking method
	public float getSpectralProbability(int score1, int score2)
	{
		float specProb = 0;

		int index1 = score1-minScore1;
		int index2 = score2-minScore2;
		
		// t1 >= score1 && t2 >= score2
		specProb += specProbWithBetterScores[index1][index2];
		
		float curCombinedSpecProb = specProb1[index1]*specProb2[index2];
		
		// t1 < score1 && t2 > score2
		int t2 = maxScore2;
		for(int t1=minScore1; t1<score1; t1++)
		{
			for(; t2>=score2; t2--)
			{
				float combinedSpecProb = specProb1[t1-minScore1]*specProb2[t2-minScore2];
				if(combinedSpecProb > curCombinedSpecProb)
					break;
			}
			if(t2+1 <= maxScore2)
				specProb += probBetterScore2[t1-minScore1][t2+1-minScore2];
		}
		
		// t1 > score1 && t2 < score2
		int t1 = maxScore1;
		for(t2=minScore2; t2<score2; t2++)
		{
			for(; t1>=score1; t1--)
			{
				float combinedSpecProb = specProb1[t1-minScore1]*specProb2[t2-minScore2];
				if(combinedSpecProb > curCombinedSpecProb)
					break;
			}
			if(t1+1 <= maxScore1)
			{
				specProb  += probBetterScore1[t1+1-minScore1][t2-minScore2];
			}
		}		
		return specProb;
	}
	
	// sum prob (t1, t2) where getNumberBetterBoth(t1,t2) <= getNumberBetterBoth(score1, score2)
	public float getNumEqualOrBetterPeptidesNumBetterPeptides(int score1, int score2)
	{
		float numBetter = 0;

		int index1 = score1-minScore1;
		int index2 = score2-minScore2;
		
		// t1 >= score1 && t2 >= score2
		float curNumRecs = numPeptidesWithBetterScores[index1][index2];
		numBetter += curNumRecs;
		
		// t1 < score1 && t2 > score2
		int t2 = maxScore2;
		for(int t1=minScore1; t1<score1; t1++)
		{
			for(; t2>=score2; t2--)
			{
				float numRecs = numPeptidesWithBetterScores[t1-minScore1][t2-minScore2];
				if(numRecs > curNumRecs)
					break;
			}
			if(t2+1 <= maxScore2)
				numBetter += numPeptidesBetterScore2[t1-minScore1][t2+1-minScore2];
		}
		
		// t1 > score1 && t2 < score2
		int t1 = maxScore1;
		for(t2=minScore2; t2<score2; t2++)
		{
			for(; t1>=score1; t1--)
			{
				float numRecs = numPeptidesWithBetterScores[t1-minScore1][t2-minScore2];
				if(numRecs > curNumRecs)
					break;
			}
			if(t1+1 <= maxScore1)
			{
				numBetter  += numPeptidesBetterScore1[t1+1-minScore1][t2-minScore2];
			}
		}		
		return numBetter;
	}	
	
	// use #peptides with t1>=score1 && t2>=score2 as the rank method
	public float getSpectralProbabilityNumBetterPeptides(int score1, int score2)
	{
		float specProb = 0;

		int index1 = score1-minScore1;
		int index2 = score2-minScore2;
		
		// t1 >= score1 && t2 >= score2
		specProb += specProbWithBetterScores[index1][index2];
		
		float curNumRecs = numPeptidesWithBetterScores[index1][index2];
		
		// t1 < score1 && t2 > score2
		int t2 = maxScore2;
		for(int t1=minScore1; t1<score1; t1++)
		{
			for(; t2>=score2; t2--)
			{
				float numRecs = numPeptidesWithBetterScores[t1-minScore1][t2-minScore2];
				if(numRecs > curNumRecs)
					break;
			}
			if(t2+1 <= maxScore2)
				specProb += probBetterScore2[t1-minScore1][t2+1-minScore2];
		}
		
		// t1 > score1 && t2 < score2
		int t1 = maxScore1;
		for(t2=minScore2; t2<score2; t2++)
		{
			for(; t1>=score1; t1--)
			{
				float numRecs = numPeptidesWithBetterScores[t1-minScore1][t2-minScore2];
				if(numRecs > curNumRecs)
					break;
			}
			if(t1+1 <= maxScore1)
			{
				specProb  += probBetterScore1[t1+1-minScore1][t2-minScore2];
			}
		}		
		return specProb;
	}
	
	public float getNumEqualOrBetterPeptidesSumScores(int score1, int score2)
	{
		float numRecs = 0;
		int sumScore = score1+score2;
		for(int s1=minScore1; s1<=maxScore1; s1++)
		{
			int s2 = sumScore-s1;
			if(s2 >= minScore2 && s2 <= maxScore2)
				numRecs += numPeptidesBetterScore1[s1-minScore1][s2-minScore2];
		}
		return numRecs;
	}
	
	public float getSpectralProbabilitySumScores(int score1, int score2)
	{
		float specProb = 0;
		int sumScore = score1+score2;
		for(int s1=minScore1; s1<=maxScore1; s1++)
		{
			int s2 = sumScore-s1;
			if(s2 >= minScore2 && s2 <= maxScore2)
			{
				specProb += probBetterScore1[s1-minScore1][s2-minScore2];
			}
		}
		return specProb;
		
	}
	
	private void initBothBetter()
	{
		numPeptidesWithBetterScores = new float[maxScore1-minScore1+1][maxScore2-minScore2+1];
		specProbWithBetterScores = new float[maxScore1-minScore1+1][maxScore2-minScore2+1];
		
		// #better scoring peptides where score1 > t1 && score2 == t2
		numPeptidesBetterScore1 = new float[maxScore1-minScore1+1][maxScore2-minScore2+1];
		probBetterScore1 = new float[maxScore1-minScore1+1][maxScore2-minScore2+1];
		
		// #better scoring peptides where score1 == t1 && score2 > t2
		numPeptidesBetterScore2 = new float[maxScore1-minScore1+1][maxScore2-minScore2+1];
		probBetterScore2 = new float[maxScore1-minScore1+1][maxScore2-minScore2+1];
		
		int maxSumScore = maxScore1+maxScore2;
		for(int energy=0; energy<=maxScore1+maxScore2-minScore1-minScore2; energy++)
		{
			for(int score1=Math.max(maxSumScore-energy-maxScore2, minScore1); score1<=maxScore1 && score1<=maxSumScore-minScore2-energy; score1++)
			{
				int score2 = maxSumScore-score1-energy;
//				System.out.println("***"+energy+" "+score1+" "+score2);
				int index1 = score1-minScore1;
				int index2 = score2-minScore2;
//				if(dist2D.getProbability(score1, score2) > 0)
//					System.out.println("=========="+score1+" "+score2+" "+" "+dist2D.getNumRecs(score1, score2)+" "+dist2D.getProbability(score1, score2));
				
				// t1 > score1 && t2 > score2
				if(score1+1 < dist2D.getMaxScore1() && score2+1 < dist2D.getMaxScore2())
				{
					numPeptidesWithBetterScores[index1][index2] += numPeptidesWithBetterScores[index1+1][index2+1];
					specProbWithBetterScores[index1][index2] += specProbWithBetterScores[index1+1][index2+1];
				}
				// t1 > score1 && t2 == score2
				if(score1+1 < dist2D.getMaxScore1())
				{
					numPeptidesWithBetterScores[index1][index2] += numPeptidesBetterScore1[index1+1][index2];
					specProbWithBetterScores[index1][index2] += probBetterScore1[index1+1][index2];
					numPeptidesBetterScore1[index1][index2] += numPeptidesBetterScore1[index1+1][index2];
					probBetterScore1[index1][index2] += probBetterScore1[index1+1][index2];
				}
				// t1 == score1 && t2 > score2
				if(score2+1 < dist2D.getMaxScore2())
				{
					numPeptidesWithBetterScores[index1][index2] += numPeptidesBetterScore2[index1][index2+1];
					specProbWithBetterScores[index1][index2] += probBetterScore2[index1][index2+1];
					numPeptidesBetterScore2[index1][index2] += numPeptidesBetterScore2[index1][index2+1];
					probBetterScore2[index1][index2] += probBetterScore2[index1][index2+1];
				}
				// t1 == score1 && t2 == score2
				float curNumRecs = dist2D.getNumRecs(score1, score2);
				numPeptidesWithBetterScores[index1][index2] += curNumRecs;
				numPeptidesBetterScore1[index1][index2] += curNumRecs;
				numPeptidesBetterScore2[index1][index2] += curNumRecs;
				
				float curProb = dist2D.getProbability(score1, score2);
				specProbWithBetterScores[index1][index2] += curProb;
				probBetterScore1[index1][index2] += curProb;
				probBetterScore2[index1][index2] += curProb;
				
//				System.out.println("specProbWithBetterScores " + score1 + " " + score2 + ": " + specProbWithBetterScores[index1][index2]);
			}
		}
		
		// init spectral probabilities
		specProb1 = new float[maxScore1-minScore1+1];
		for(int t = maxScore1; t>=minScore1; t--)
		{
			specProb1[t-minScore1] = probBetterScore2[t-minScore1][0];
			if(t+1 <= maxScore1)
				specProb1[t-minScore1] += specProb1[t+1-minScore1];
		}
		specProb2 = new float[maxScore2-minScore2+1];
		for(int t = maxScore2; t>=minScore2; t--)
		{
			specProb2[t-minScore2] = probBetterScore1[0][t-minScore2];
			if(t+1 <= maxScore2)
				specProb2[t-minScore2] += specProb2[t+1-minScore2];
		}
		
		/*
		float sum = 0;
		for(int s1=minScore1; s1<=maxScore1; s1++)
		{
			for(int s2=minScore2; s2<=maxScore2; s2++)
			{
//				if(s1 < 110 && s2 < 90)
//					continue;
//				if(dist2D.getNumRecs(s1, s2) > 0)
				if(s1+s2 >= 56-3)
				{
//					if(numPeptidesWithBetterScores[s1-minScore1][s2-minScore2] <= 4)
						sum += dist2D.getProbability(s1, s2);
					System.out.println(s1+" "+s2+" "+dist2D.getNumRecs(s1, s2)+" "+numPeptidesWithBetterScores[s1-minScore1][s2-minScore2]
							+" "+dist2D.getProbability(s1, s2)+" "+specProbWithBetterScores[s1-minScore1][s2-minScore2]);
				}
			}
		}
		System.out.println("Sum: "+sum);
		System.exit(0);
		*/
		
	}
	
	// If one wants to compute spec prob. of all score pairs, call the following method
	// specProb[i][j]: sum of probabilities of score pairs (t1, t2) where numBetterPeptides of (t1, t2) <= numBetterPeptides of (i,j)
	private float[][] specProbMerged;
	// numBetterPeptidesMerged[i][j]: number of peptides at score pairs (t1, t2) where numBetterPeptides of (t1, t2) <= numBetterPeptides of (i,j)
	private float[][] numBetterPeptidesMerged;
	protected void initMerged()
	{
		numBetterPeptidesMerged = new float[maxScore1-minScore1+1][maxScore2-minScore2+1];
		specProbMerged = new float[maxScore1-minScore1+1][maxScore2-minScore2+1];
		
		// initialization
		numBetterPeptidesMerged[dist2D.getMaxScore1()-1-minScore1][dist2D.getMaxScore2()-1-minScore2] = dist2D.getNumRecs(dist2D.getMaxScore1()-1,dist2D.getMaxScore2()-1);
		specProbMerged[dist2D.getMaxScore1()-1-minScore1][dist2D.getMaxScore2()-1-minScore2] = dist2D.getProbability(dist2D.getMaxScore1()-1,dist2D.getMaxScore2()-1);
		
		int minScore1 = dist2D.getMinScore1();
		int minScore2 = dist2D.getMinScore2();
		int maxScore1 = dist2D.getMaxScore1()-1;
		int maxScore2 = dist2D.getMaxScore2()-1;
		
		int maxSumScore = dist2D.getMaxScore1()+dist2D.getMaxScore2()-2;
		for(int energy=1; energy<dist2D.getMaxScore1()+dist2D.getMaxScore2()-1; energy++)
		{
			for(int score1=Math.max(maxSumScore-energy-maxScore2, minScore1); score1<=maxScore1 && score1<=maxSumScore-minScore2-energy; score1++)
			{
				int score2 = maxSumScore-score1-energy;
				
				int index1 = score1-minScore1;
				int index2 = score2-minScore2;
				
				// t1 >= score1 && t2 >= score2
				numBetterPeptidesMerged[index1][index2] += numPeptidesWithBetterScores[index1][index2];
				specProbMerged[index1][index2] += specProbWithBetterScores[index1][index2];
				
				float curNumRecs = numPeptidesWithBetterScores[index1][index2];
				
				// t1 < score1 && t2 > score2
				int t2 = maxScore2;
				for(int t1=minScore1; t1<score1; t1++)
				{
					for(; t2>=score2; t2--)
					{
						float numRecs = numPeptidesWithBetterScores[t1-minScore1][t2-minScore2];
						if(numRecs > curNumRecs)
							break;
					}
					if(t2+1 <= maxScore2)
					{
						numBetterPeptidesMerged[index1][index2] += numPeptidesBetterScore2[t1-minScore1][t2+1-minScore2];
						specProbMerged[index1][index2] += probBetterScore1[t1-minScore1][t2+1-minScore2];
					}
				}
				
				// t1 > score1 && t2 < score2
				int t1 = maxScore1;
				for(t2=minScore2; t2<score2; t2++)
				{
					for(; t1>=score1; t1--)
					{
						float numRecs = numPeptidesWithBetterScores[t1-minScore1][t2-minScore2];
						if(numRecs > curNumRecs)
							break;
					}
					if(t1+1 <= maxScore1)
					{
						numBetterPeptidesMerged[index1][index2] += numPeptidesBetterScore1[t1+1-minScore1][t2-minScore2];
						specProbMerged[index1][index2] += probBetterScore1[t1+1-minScore1][t2-minScore2];
					}
				}
				// t1 < score1 && t2 < score2: don't need to be considered
			}
		}
	}
	
}
