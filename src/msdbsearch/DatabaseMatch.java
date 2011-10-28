package msdbsearch;

import java.util.Comparator;

import msgf.ScoreDist;

public class DatabaseMatch implements Comparable<DatabaseMatch> {
	private int		index;
	private byte 	length;
	private int 	score;
	private int		nominalPeptideMass;
	
	// optional
	private int		deNovoScore;	
	private double	specProb;
	private ScoreDist	scoreDist;
	private boolean isProteinNTerm;
	private boolean isProteinCTerm;
	private int		charge;

	private String pepSeq;
	public DatabaseMatch(int index, int length, int score, int nominalPeptideMass, String pepSeq) {
		super();
		this.index = index;
		this.length = (byte)length;
		this.score = score;
		this.nominalPeptideMass = nominalPeptideMass;
		this.pepSeq = pepSeq;
		isProteinNTerm = false;
		isProteinCTerm = false;
	}
	
	public DatabaseMatch setProteinNTerm(boolean isProteinNTerm)
	{
		this.isProteinNTerm = isProteinNTerm;
		return this;
	}

	public DatabaseMatch setProteinCTerm(boolean isProteinCTerm)
	{
		this.isProteinCTerm = isProteinCTerm;
		return this;
	}

	public DatabaseMatch setCharge(int charge)
	{
		this.charge = charge;
		return this;
	}
	
	public int getNominalPeptideMass()
	{
		return nominalPeptideMass;
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
	
	public int getCharge() {
		return charge;
	}

	public boolean isProteinNTerm()
	{
		return isProteinNTerm;
	}

	public boolean isProteinCTerm()
	{
		return isProteinCTerm;
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

	public void setSpecProb(double specProb)
	{
		this.specProb = specProb;
	}
	
	public double getSpecProb() {
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
	
	public static class SpecProbComparator implements Comparator<DatabaseMatch>
	{
		@Override
		public int compare(DatabaseMatch arg0, DatabaseMatch arg1) {
			if(arg0.getSpecProb() < arg1.getSpecProb())
				return 1;
			else if(arg0.getSpecProb() > arg1.getSpecProb())
				return -1;
			else
				return 0;
		}
		
	}
	
}
