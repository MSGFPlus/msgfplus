package misc;

import java.io.BufferedReader;
import java.io.FileReader;

import msgf.Histogram;
import msutil.AminoAcidSet;

public class CalcFastaDBSize {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 1 || (!argv[0].endsWith(".fasta") && !argv[0].endsWith(".fa")))
		{
			System.out.println("usage: java CalcFastaDBSize *.fasta");
			System.exit(0);
		}
		BufferedReader in = new BufferedReader(new FileReader(argv[0]));
		String s;
		int length = 0;
		int numProteins = 0;
		int[] numTrypticPeptides = new int[11];	// the number of tryptic peptides allowing no miscleavages
		int pepLen = 0;
		Histogram<Character> aaHist = new Histogram<Character>();
//		StringBuffer pep = new StringBuffer();
//		java.util.HashSet<String> pepSet = new java.util.HashSet<String>();
		while((s=in.readLine()) != null)
		{
			if(s.startsWith(">"))
			{
				numProteins++;
				pepLen = 0;
//				pep = new StringBuffer();
				continue;
			}
			else
			{
				length += s.length();
				for(int i=0; i<s.length(); i++)
				{
					char c = s.charAt(i);
					if(Character.isLetter(c))
						aaHist.add(c);
					if(c == 'K' || c == 'R')
					{
						if(pepLen >= 4 && pepLen < numTrypticPeptides.length-1)
						{
//							pep.append(c);
//							if(pep.length() == 5)
//								System.out.println(pep);
							numTrypticPeptides[++pepLen]++;
//							if(pepLen == 10)
//								pepSet.add(pep.toString());
						}
						pepLen = 0;
//						pep = new StringBuffer();
					}
					else if(AminoAcidSet.getStandardAminoAcidSet().getAminoAcid(c) == null)
					{
						pepLen = 0;
//						pep = new StringBuffer();
					}
					else
					{
						pepLen++;
//						pep.append(c);
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
//		System.out.println("Distinct Peptides length 10: " + pepSet.size());
		System.out.println("Amino Acid Composition:");
		aaHist.printSorted();
	}
}
