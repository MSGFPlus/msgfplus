package edu.ucsd.msjava.misc;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Peptide;

public class MassCalc {
	public static void main(String argv[])
	{
		AminoAcidSet aaSet;
		if(argv.length > 0)
			aaSet = AminoAcidSet.getStandardAminoAcidSet();
		else
			aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		while(true)
		{
			System.out.print("Sequence : ");
			String pepStr = readLine();
			if(pepStr.matches("[A-Z]*\\..+\\.[A-Z]*"))
				pepStr = pepStr.substring(pepStr.indexOf('.')+1, pepStr.lastIndexOf('.'));

			Peptide pep = new Peptide(pepStr, aaSet);
			Composition composition = null;
			if(!pep.isModified())
				composition = pep.getComposition();
			System.out.println("Mass : " + pep.getAccurateMass() + " (" + composition + ")\t" 
					+ pep.getNominalMass() + "\t"
					+ (pep.getMass() + Composition.H2O) + "\t"
					+ (pep.getMass() + Composition.H2O + Composition.PROTON) + "\t"
					+ ((pep.getMass() + Composition.H2O + 2*Composition.PROTON)/2) + "\t" 
					+ ((pep.getMass() + Composition.H2O + 3*Composition.PROTON)/3) 
					);
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