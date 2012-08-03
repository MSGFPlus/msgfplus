/***************************************************************************
  * Title:          
  * Author:         Sangtae Kim
  * Last modified:  
  *
  * Copyright (c) 2008-2009 The Regents of the University of California
  * All Rights Reserved
  * See file LICENSE for details.
  ***************************************************************************/
package edu.ucsd.msjava.fdr;

/**
 * The Class ScoredString.
 */
public class ScoredString extends Pair<String, Float> implements Comparable<Pair<String, Float>>
{
	
	/**
	 * Instantiates a new scored string.
	 * 
	 * @param peptide the peptide
	 * @param score the score
	 */
	public ScoredString(String peptide, Float score)
	{
		super(peptide, score);
	}
	
	/**
	 * Instantiates a new scored string, using an integer score.
	 * @param score
	 * @param peptide
	 */
	public ScoredString(String peptide, int score)
	{
		super(peptide, (float)score);
	}
	
	public int compareTo(Pair<String, Float> o) {
		int scoreComp = getSecond().compareTo(o.getSecond());
		if(scoreComp != 0)
			return scoreComp;
		else
			return getFirst().compareTo(o.getFirst());
	}
	
	/**
	 * Gets the str.
	 * 
	 * @return the str
	 */
	public String getStr() {
		return super.getFirst();
	}

	/**
	 * Gets the score.
	 * 
	 * @return the score
	 */
	public float getScore() {
		return super.getSecond();
	}
	
}

