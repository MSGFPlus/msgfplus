package msdictionary;

import java.util.Hashtable;

public class Codon {
	private static Hashtable<String, Character> table;
	
	public static char translate(String codon)
	{
		Character aa = table.get(codon);
		if(aa == null)
			return '\0';
		else
			return aa;
	}

	public static char translateRevComplement(String codon)
	{
		StringBuffer revComp = new StringBuffer();
		for(int i=codon.length()-1; i>=0; i--)
			revComp.append(complement(codon.charAt(i)));
		Character aa = table.get(revComp.toString());
		if(aa == null)
			return '\0';
		else
			return aa;
	}
	
	public static char complement(char aa)
	{
		if(aa == 'A')
			return 'T';
		else if(aa == 'T')
			return 'A';
		else if(aa == 'G')
			return 'C';
		else if(aa == 'C')
			return 'G';
		else
			return '\0';
	}
	static {
		table = new Hashtable<String, Character>();
		table.put("ATT", 'I');
		table.put("ATC", 'I');
		table.put("ATA", 'I');
		
		table.put("CTT", 'L');
		table.put("CTC", 'L');
		table.put("CTA", 'L');
		table.put("CTG", 'L');
		table.put("TTA", 'L');
		table.put("TTG", 'L');

		table.put("GTT", 'V');
		table.put("GTC", 'V');
		table.put("GTA", 'V');
		table.put("GTG", 'V');
		
		table.put("TTT", 'F');
		table.put("TTC", 'F');

		table.put("ATG", 'M');
		
		table.put("TGT", 'C');
		table.put("TGC", 'C');

		table.put("GCT", 'A');
		table.put("GCC", 'A');
		table.put("GCA", 'A');
		table.put("GCG", 'A');

		table.put("GGT", 'G');
		table.put("GGC", 'G');
		table.put("GGA", 'G');
		table.put("GGG", 'G');

		table.put("CCT", 'P');
		table.put("CCC", 'P');
		table.put("CCA", 'P');
		table.put("CCG", 'P');
		
		table.put("ACT", 'T');
		table.put("ACC", 'T');
		table.put("ACA", 'T');
		table.put("ACG", 'T');

		table.put("TCT", 'S');
		table.put("TCC", 'S');
		table.put("TCA", 'S');
		table.put("TCG", 'S');
		table.put("AGT", 'S');
		table.put("AGC", 'S');

		table.put("TAT", 'Y');
		table.put("TAC", 'Y');

		table.put("TGG", 'W');
		
		table.put("CAA", 'Q');
		table.put("CAG", 'Q');

		table.put("AAT", 'N');
		table.put("AAC", 'N');

		table.put("CAT", 'H');
		table.put("CAC", 'H');

		table.put("GAA", 'E');
		table.put("GAG", 'E');
		
		table.put("GAT", 'D');
		table.put("GAC", 'D');

		table.put("AAA", 'K');
		table.put("AAG", 'K');

		table.put("CGT", 'R');
		table.put("CGC", 'R');
		table.put("CGA", 'R');
		table.put("CGG", 'R');
		table.put("AGA", 'R');
		table.put("AGG", 'R');

		table.put("TAA", 'X');
		table.put("TAG", 'X');
		table.put("TGA", 'X');
	}
}
