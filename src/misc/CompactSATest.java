package misc;

import msdbsearch.CompactFastaSequence;
import msdbsearch.CompactSuffixArray;

public class CompactSATest {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 1)
			printUsageAndExit();
		readTest(argv[0]);
	}
	
	public static void printUsageAndExit()
	{
		System.out.println("Usage: java CompactSATest *.fasta (or *.fa)");
		System.exit(-1);
	}
	
	public static void readTest(String fileName) throws Exception
	{
//		String fileName = "/home/sangtaekim/Research/Data/RNASeq/AminPark_tryptic_pep"+DECOY_DB_EXTENSION;
		CompactSuffixArray sa = new CompactSuffixArray(new CompactFastaSequence(fileName), 30);
		System.out.println("Read complete");
	}
}
