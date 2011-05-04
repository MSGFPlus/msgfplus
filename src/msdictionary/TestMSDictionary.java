package msdictionary;

import suffixarray.SuffixArraySequence;
import suffixarray.SuffixArray;

public class TestMSDictionary {
	public static void main(String argv[]) throws Exception
	{
//		testSuffixArrayConstruction();
//		translateAndBuildSuffixArray();
	}
	
	public static void testDictionary() throws Exception
	{
	}
	
	public static void testSuffixArrayConstruction()
	{
		int splitNum=0;
		String fileName = System.getProperty("user.home")+"/Research/Data/HumanGenome/translated/HSRM.NCBI36.54.translation."+splitNum+".fasta";
		SuffixArraySequence sequence = new SuffixArraySequence(fileName);
		System.out.println("FastaSequence done");
		SuffixArray sa = new SuffixArray(sequence);
		System.out.println("SuffixArray done");
		System.out.println(sa.search("PDLHSSYALARAWAGHGSHERAQS"));
		System.out.println(sa.search("KPVRKGIPANARLPMVRE"));
		System.out.println(sa.search("DMFAVSMMWSLEAKPVRKGIPANARLPMVRE"));
		/*
//		String fileName = System.getProperty("user.home")+"/Research/Data/SProt/test2/uniprot_sprot.fasta";
		for(int splitNum=2; splitNum<4; splitNum++)
		{
			String fileName = System.getProperty("user.home")+"/Research/Data/HumanGenome/translated/HSRM.NCBI36.54.translation."+splitNum+".fasta";
			Adapter sequence = new FastaSequence(fileName);
			System.out.println("FastaSequence done");
			SuffixArray sa = new SuffixArray(sequence);
			System.out.println("SuffixArray done");
		}
		*/
	}
	
	public static void translateAndBuildSuffixArray()
	{
		for(int splitNum=0; splitNum<4; splitNum++)
		{
			String genomeFileName = System.getProperty("user.home")+"/Research/Data/HumanGenome/splitted/Homo_sapiens.NCBI36.54.dna_rm."+splitNum+".fasta";
			String translationFileName = System.getProperty("user.home")+"/Research/Data/HumanGenome/translated/HSRM.NCBI36.54.translation."+splitNum+".fasta";
			
			// translation
			new GenomeTranslator(genomeFileName).translateAndWriteTo(translationFileName);
			
			// Build Suffix tree
			SuffixArraySequence sequence = new SuffixArraySequence(translationFileName);
			System.out.println("FastaSequence done");
			SuffixArray sa = new SuffixArray(sequence);
			System.out.println("SuffixArray done " + sa);
		}
		
	}
}
