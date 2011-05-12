package msdbsearch;

import msgf.ScoreDist;

public class DatabaseMatch implements Comparable<DatabaseMatch> {
	private int		index;
	private byte 	length;
	private int 	score;
	
	// optional
	private int		deNovoScore;	
	private float	specProb;
	private ScoreDist	scoreDist;
	private boolean isProteinNTerm;

	private String pepSeq;
	public DatabaseMatch(int index, int length, int score) {
		super();
		this.index = index;
		this.length = (byte)length;
		this.score = score;
		pepSeq = null;
		isProteinNTerm = false;
	}
	
	public DatabaseMatch pepSeq(String pepSeq)
	{
		this.pepSeq = pepSeq;
		return this;
	}

	public DatabaseMatch setProteinNTerm(boolean isProteinNTerm)
	{
		this.isProteinNTerm = isProteinNTerm;
		return this;
	}
	
	public String getPepSeq()
	{
		return pepSeq;
	}
	
	public int getIndex() {
		return index;
	}

	public int getLength() {
		return length;
	}

	public boolean isProteinNTerm()
	{
		return isProteinNTerm;
	}
	
	public void setScore(int score)
	{
		this.score = score;
	}
	
	public int getScore() {
		return score;
	}
	
	public void setDeNovoScore(int deNovoScore)
	{
		this.deNovoScore = deNovoScore;
	}
	
	public int getDeNovoScore()	{
		return deNovoScore;
	}

	public void setSpecProb(float specProb)
	{
		this.specProb = specProb;
	}
	
	public float getSpecProb() {
		return specProb;
	}
	
	public void setScoreDist(ScoreDist scoreDist)
	{
		this.scoreDist = scoreDist;
	}
	
	public ScoreDist getScoreDist()
	{
		return scoreDist;
	}
	
	public int hashCode()	{
		return index*length;
	}
	
	public boolean equals(Object obj) {
		if(obj instanceof DatabaseMatch)
		{
			DatabaseMatch other = (DatabaseMatch)obj;
			if(index == other.index && length == other.length)
				return true;
		}
		return false;
	}
	
	@Override
	public int compareTo(DatabaseMatch o) {
		float diff = score - o.score;
		if(diff > 0)
			return 1;
		else if(diff == 0)
		{
			int diff2 = index-o.index;
			if(diff2 != 0)
				return diff2;
			else
				return length-o.length;
		}
		else
			return -1;
	}
	
}
