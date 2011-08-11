package misc;

import suffixarray.SuffixArraySequence;
import msdbsearch.SuffixArrayForMSGFDB;

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
//		String fileName = "/home/sangtaekim/Research/Data/RNASeq/AminPark_tryptic_pep.revConcat.fasta";
		SuffixArrayForMSGFDB sa = new SuffixArrayForMSGFDB(new SuffixArraySequence(fileName), 6, 30);
		System.out.println("Read complete");
	}
}
