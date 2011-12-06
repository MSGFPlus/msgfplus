package misc;

import msdbsearch.CompactFastaSequence;
import msdbsearch.CompactSuffixArray;
import msutil.AminoAcid;
import msutil.AminoAcidSet;

public class MSGFPlusPaper {
	public static void main(String argv[]) throws Exception
	{
//		nominalMassTable();
		checkPeptidesWithNominalMassErrors();
	}
	
	public static void nominalMassTable() throws Exception
	{
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		aaSet = AminoAcidSet.getStandardAminoAcidSet();
		for(AminoAcid aa : aaSet)
		{
			float mass = aa.getMass();
			int nominalMass = aa.getNominalMass();
			float rescaledMass = mass*0.9995f;
			float error = rescaledMass-nominalMass;
			float errorPPM = (rescaledMass-nominalMass)/mass*1e6f;
			System.out.println(aa.getResidue()+"\t"+mass+"\t"+nominalMass+"\t"+rescaledMass+"\t"+error+"\t"+errorPPM);
		}
	}
	
	public static void checkPeptidesWithNominalMassErrors() throws Exception
	{
		String fileName = "/home/sangtaekim/Research/Data/IPI/IPI_human_3.87.fasta";
		int maxPeptideLength = 40;
		CompactFastaSequence fastaSequence = new CompactFastaSequence(fileName);
		CompactSuffixArray sa = new CompactSuffixArray(fastaSequence, maxPeptideLength);
		sa.measureNominalMassError(AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys());
	}
}
