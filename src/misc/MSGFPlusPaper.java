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
		String delimiter = " & ";
		System.out.println("Residue & Mass & NominalMass & RescaledMass & RescalingError & RescalingErrPPM\\\\");
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
//		aaSet = AminoAcidSet.getStandardAminoAcidSet();
		float sum = 0;
		float sumPPM = 0;
		for(AminoAcid aa : aaSet)
		{
			float mass = aa.getMass();
			int nominalMass = aa.getNominalMass();
			float rescaledMass = mass*0.9995f;
			float error = rescaledMass-nominalMass;
			float errorPPM = (rescaledMass-nominalMass)/mass*1e6f;
//			System.out.println(aa.getResidue()+delimiter+mass+delimiter+nominalMass+delimiter+rescaledMass+delimiter+error+delimiter+errorPPM+"\\\\");
			System.out.format("%c%s%.3f%s%d%s%.3f%s%.3f%s%.3f\\\\\n",
					aa.getResidue(),delimiter,
					mass,delimiter,
					nominalMass,delimiter,
					rescaledMass,delimiter,
					error,delimiter,
					errorPPM);
			sum += error;
			sumPPM += errorPPM;
		}
		System.out.println("AverageError\t"+sum/20);
		System.out.println("AverageErrorPPM\t"+sumPPM/20);
	}
	
	public static void checkPeptidesWithNominalMassErrors() throws Exception
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/IPI/IPI_human_3.87.fasta";
		int maxPeptideLength = 100;
		CompactFastaSequence fastaSequence = new CompactFastaSequence(fileName);
		CompactSuffixArray sa = new CompactSuffixArray(fastaSequence, maxPeptideLength);
		sa.measureNominalMassError(AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys());
	}
}
