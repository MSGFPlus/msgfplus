package msutil;

public class Annotation {
	private AminoAcid prevAA;
	private Peptide peptide;
	private AminoAcid nextAA;
	
	public Annotation(AminoAcid prevAA, Peptide peptide, AminoAcid nextAA)
	{
		this.prevAA = prevAA;
		this.peptide = peptide;
		this.nextAA = nextAA;
	}
	
	public Annotation(String annotationStr, AminoAcidSet aaSet)
	{
		String pepStr = annotationStr.substring(annotationStr.indexOf('.')+1, annotationStr.lastIndexOf('.'));
		char prevAAResidue = annotationStr.charAt(0);
		char nextAAResidue = annotationStr.charAt(annotationStr.length()-1);
		
		prevAA = aaSet.getAminoAcid(prevAAResidue);
		peptide = aaSet.getPeptide(pepStr);
		nextAA = aaSet.getAminoAcid(nextAAResidue);
	}
	
	public boolean isProteinNTerm()
	{
		return prevAA == null;
	}
	
	public boolean isProteiNCTerm()
	{
		return nextAA == null;
	}
	
	public AminoAcid getPrevAA() {
		return prevAA;
	}

	public void setPrevAA(AminoAcid prevAA) {
		this.prevAA = prevAA;
	}

	public Peptide getPeptide() {
		return peptide;
	}

	public void setPeptide(Peptide peptide) {
		this.peptide = peptide;
	}

	public AminoAcid getNextAA() {
		return nextAA;
	}

	public void setNextAA(AminoAcid nextAA) {
		this.nextAA = nextAA;
	}
	
	@Override
	public String toString()
	{
		if(peptide == null)
			return null;
		StringBuffer output = new StringBuffer();
		if(prevAA != null)
			output.append(prevAA.getResidueStr());
		output.append("."+peptide.toString()+".");
		if(nextAA != null)
			output.append(nextAA.getResidueStr());
		return output.toString();
	}
}
