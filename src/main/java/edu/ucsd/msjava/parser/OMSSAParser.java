package edu.ucsd.msjava.parser;

import java.io.FileNotFoundException;
import java.util.ArrayList;

import edu.ucsd.msjava.msutil.Peptide;


public class OMSSAParser {
	public static PSMList<PSM> parse(String fileName)
	{
		BufferedLineReader in = null;
		try {
			in = new BufferedLineReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		PSMList<PSM> matches = new PSMList<PSM>();
		
		String s;
		in.readLine();	// ignore the first line
		while((s=in.readLine()) != null)
		{
			String[] token = splitCSVLine(s);
			if(token.length != 15)
				continue;
			
			int specNum = Integer.parseInt(token[0]);
			String title = token[1];
			Peptide peptide = new Peptide(token[2]);	// amino acid set
			float eValue = Float.parseFloat(token[3]);
			String protein = token[9];
			String mod = token[10];
			int charge = Integer.parseInt(token[11]);
			
			matches.add(new PSM().scanNum(specNum).title(title).peptide(peptide)
					.probScore(eValue).protein(protein).ptm(mod).charge(charge));
		}
		return matches;
	}
	
	private static String[] splitCSVLine(String s)
	{
		ArrayList<String> token = new ArrayList<String>();
		StringBuffer buf = new StringBuffer();
		boolean ignoreComma = false;
		for(int i=0; i<s.length(); i++)
		{
			char c = s.charAt(i);
			if(!ignoreComma && c == ',')
			{
				token.add(buf.toString().trim());
				buf = new StringBuffer();
			}
			else if(c == '"')
			{
				ignoreComma = !ignoreComma;
			}
			else
				buf.append(c);
		}
		token.add(buf.toString().trim());
		return token.toArray(new String[0]);
	}
	
	public static void main(String argv[])
	{
		test();
	}
	
	public static void test()
	{
		String[] targetResult = {
				System.getProperty("user.home")+"/Research/Data/Heck/omssa_ETD_090309_sm4067_01_sprot.csv",
				System.getProperty("user.home")+"/Research/Data/Heck/omssa_ETD_090309_sm4067_03_sprot.csv",
		};
		String[] decoyResult = {
				System.getProperty("user.home")+"/Research/Data/Heck/omssa_ETD_090309_sm4067_01_revsprot.csv",
				System.getProperty("user.home")+"/Research/Data/Heck/omssa_ETD_090309_sm4067_03_revsprot.csv",
		};
		
		PSMList<PSM> targetPSM = OMSSAParser.parse(targetResult[1]).getDistinctiveSpectralSet().getDistinctivePeptideSet();
		PSMList<PSM> decoyPSM = OMSSAParser.parse(decoyResult[1]).getDistinctiveSpectralSet().getDistinctivePeptideSet();
		PSMList<PSM> significantPSMList = PSMList.selectUsingFDR(targetPSM, decoyPSM, 0.05f);
		@SuppressWarnings("unused")
    PSMList<PSM> significantPSMPeptideList = significantPSMList.getDistinctivePeptideSet();
		int size = 0;
		for(PSM psm : significantPSMList)
		{
			System.out.println(psm);
			size++;
		}
		System.out.println("Size: " + size);
		
	}
}
