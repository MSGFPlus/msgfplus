package edu.ucsd.msjava.msgf;

import java.util.HashMap;
import java.util.Map.Entry;

import edu.ucsd.msjava.msutil.Annotation;
import edu.ucsd.msjava.msutil.Matter;

public class GeneratingFunctionGroup<T extends Matter> extends HashMap<T, GeneratingFunction<T>> implements GF<T> {

	private static ScoreDistFactory factory = new ScoreDistFactory(false, true);
	private static final long serialVersionUID = 1L;
	private ScoreDist mergedScoreDist = null;
	
	public void registerGF(T sink, GeneratingFunction<T> gf)
	{
		this.put(sink, gf);
	}
	
	public boolean computeGeneratingFunction()
	{
		int minScore = Integer.MAX_VALUE;
		int maxScore = Integer.MIN_VALUE;
		for(Entry<T, GeneratingFunction<T>> entry : this.entrySet())
		{
			GeneratingFunction<T> gf = entry.getValue();
			if(!gf.isGFComputed())
			{
				if(gf.computeGeneratingFunction() == true)
				{
					int curMinScore = gf.getMinScore();
					if(minScore > curMinScore)
						minScore = curMinScore;
					int curMaxScore = gf.getMaxScore();
					if(maxScore < curMaxScore)
						maxScore = curMaxScore;
				}
			}
		}
		if(minScore >= maxScore)
			return false;
		mergedScoreDist = factory.getInstance(minScore, maxScore);
		for(Entry<T, GeneratingFunction<T>> entry : this.entrySet())
		{
			GeneratingFunction<T> gf = entry.getValue();
			mergedScoreDist.addProbDist(gf.getScoreDist(), 0, 1f);
		}		
		return true;
	}
	
	public int getScore(Annotation annotation)
	{
		int score = Integer.MIN_VALUE;
		for(Entry<T, GeneratingFunction<T>> entry : this.entrySet())
		{
			GeneratingFunction<T> gf = entry.getValue();
			int curScore = gf.getScore(annotation);
			if(curScore > score)
				score = curScore;
		}		
		
		return score;
	}
	
	public double getSpectralProbability(int score)
	{
		return mergedScoreDist.getSpectralProbability(score);
	}
	
	public int getMaxScore()
	{
		if(mergedScoreDist == null)
			System.out.println("Debug");
		return mergedScoreDist.getMaxScore();
	}
	
	public ScoreDist getScoreDist()
	{
		return mergedScoreDist;
	}
}
