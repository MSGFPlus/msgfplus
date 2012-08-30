package edu.ucsd.msjava.msdbsearch;

import java.util.TreeSet;
import java.util.List;
import java.util.SortedSet;

import edu.ucsd.msjava.msutil.ActivationMethod;

public class DatabaseMatch extends Match {
	private int		index;
	private byte 	length;
	
	// optional
	private boolean isProteinNTerm;
	private boolean isProteinCTerm;
	private boolean isNTermMetCleaved = false;
	
	private Float psmQValue = null;
	private Float pepQValue = null;
	
	// for degenerate peptides
	private SortedSet<Integer> indices;
	
	public DatabaseMatch(
			int index, 
			byte length,
			int score, 
			float peptideMass, 
			int nominalPeptideMass,
			int charge, 
			String pepSeq,
			ActivationMethod[] actMethodArr
			) 
	{
		super(score, peptideMass, nominalPeptideMass, charge, pepSeq, actMethodArr);
		this.index = index;
		this.length = length;
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

	public DatabaseMatch setNTermMetCleaved(boolean isNTermMetCleaved)
	{
		this.isNTermMetCleaved = isNTermMetCleaved;
		return this;
	}
	
	public boolean isNTermMetCleaved()
	{
		return this.isNTermMetCleaved;
	}
	
	public void setPSMQValue(float psmQValue)
	{
		this.psmQValue = psmQValue;
	}

	public Float getPSMQValue()	
	{
		return this.psmQValue;
	}
	
	public void setPepQValue(Float pepQValue)
	{
		this.pepQValue = pepQValue;
	}
	
	public Float getPepQValue()
	{
		return this.pepQValue;
	}
	public void addIndex(int index)
	{
		if(indices == null)
		{
			indices = new TreeSet<Integer>();
			indices.add(this.index);
		}
		indices.add(index);
	}
	
	public SortedSet<Integer> getIndices()
	{
		if(indices == null)
		{
			SortedSet<Integer> temp = new TreeSet<Integer>();
			temp.add(index);
			return temp;
		}
		return indices;
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

	public boolean isProteinCTerm()
	{
		return isProteinCTerm;
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
}
