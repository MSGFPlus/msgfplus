package misc;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.Peptide;

public class MassCalc {
	public static void main(String argv[])
	{
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		if(argv.length == 0)
		{
			while(true)
			{
				System.out.print("Sequence : ");
				Peptide seq = aaSet.getPeptide(readLine());
				Composition composition = seq.getComposition();
				System.out.println("Mass : " + seq.getAccurateMass() + " (" + composition + ")\t" 
						+ (seq.getMass() + Composition.H2O) + "\t"
						+ (seq.getMass() + Composition.H2O + Composition.PROTON) + "\t"
						+ ((seq.getMass() + Composition.H2O + 2*Composition.PROTON)/2) + "\t" 
						+ ((seq.getMass() + Composition.H2O + 3*Composition.PROTON)/3) 
						);
			}
		}
		else
		{
			if(argv[0].equalsIgnoreCase("0"))
			{
				while(true)
				{
					System.out.print("Sequence : ");
					Peptide seq = aaSet.getPeptide(readLine());
					float mass = 0;
					for(AminoAcid aa : seq)
					{
						mass += aa.getNominalMass();
					}
					System.out.println("Mass : " + mass);
				}
			}
		}
	}
	
	public static String readLine()
	{
		String buffer = null;
		
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
			buffer = br.readLine();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return buffer;
	}
	
}