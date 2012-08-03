package edu.ucsd.msjava.scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;

import edu.ucsd.msjava.msgf.Histogram;
import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;

public class GetDBInfo {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 1 || (!argv[0].endsWith(".fasta") && !argv[0].endsWith(".fa")))
		{
			System.err.println("usage: java CalcFastaDBSize *.fasta");
			System.exit(-1);
		}
		File dbFile = new File(argv[0]);
		if(!dbFile.exists())
		{
			System.err.println(argv[0] + " doen't exist!");
			System.exit(-1);
		}
		getDBInfo(dbFile);
	}
	
	private static final int MAX_PEPTIDE_LENGTH = 50;
	
	public static void getDBInfo(File dbFile) throws Exception
	{
		BufferedReader in = new BufferedReader(new FileReader(dbFile));
		String s;
		int length = 0;
		int numProteins = 0;
		int[] numTrypticPeptides = new int[MAX_PEPTIDE_LENGTH+1];	// the number of tryptic peptides allowing no miscleavages
		Histogram<Integer> hist = new Histogram<Integer>();
		
		int pepLen = 0;
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		while((s=in.readLine()) != null)
		{
			if(s.startsWith(">"))
			{
				numProteins++;
				pepLen = 0;
				continue;
			}
			else
			{
				length += s.length();
				for(int i=0; i<s.length(); i++)
				{
					char c = s.charAt(i);
					AminoAcid aa = aaSet.getAminoAcid(c);
					if(aa != null)
						hist.add(aaSet.getIndex(aa));
					if(c == 'K' || c == 'R')
					{
						if(pepLen >= 4 && pepLen < numTrypticPeptides.length-1)
						{
							numTrypticPeptides[++pepLen]++;
						}
						pepLen = 0;
					}
					else if(aa == null)
					{
						pepLen = 0;
					}
					else
					{
						pepLen++;
					}
				}
			}
		}
		System.out.println("#Proteins: " + numProteins);
		System.out.println("#Amino acids: " + length);
		System.out.println("#Tryptic Peptides (with no miscleavage)");
		System.out.println("Length\tNumber");
		for(int i=5; i<numTrypticPeptides.length; i++)
			System.out.println(i+"\t"+numTrypticPeptides[i]);
		System.out.println("Amino acid composition:");
		ArrayList<Integer> keyList = new ArrayList<Integer>(hist.keySet());
		Collections.sort(keyList);
		for(Integer key : keyList)
			System.out.println(aaSet.getAminoAcid(key).getResidueStr()+"\t"+(hist.get(key)/(float)hist.totalCount()));
	}
}
