package msdbsearch;

public class DatabaseMatch extends Match {
	private int		index;
	private byte 	length;
	
	// optional
	private boolean isProteinNTerm;
	private boolean isProteinCTerm;

	public DatabaseMatch(
			int index, 
			byte length,
			int score, 
			float peptideMass, 
			int nominalPeptideMass,
			int charge, 
			String pepSeq 
			) 
	{
		super(score, peptideMass, nominalPeptideMass, charge, pepSeq);
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
