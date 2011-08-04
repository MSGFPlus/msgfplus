package misc;

import suffixarray.SuffixArraySequence;
import msdbsearch.SuffixArrayForMSGFDB;

public class CompactSATest {
	public static void main(String argv[]) throws Exception
	{
		readTest();
	}
	
	public static void readTest() throws Exception
	{
		String fileName = "/home/sangtaekim/Research/Data/RNASeq/AminPark_tryptic_pep.revConcat.fasta";
		SuffixArrayForMSGFDB sa = new SuffixArrayForMSGFDB(new SuffixArraySequence(fileName), 6, 30);
		System.out.println("Read complete");
	}
}
