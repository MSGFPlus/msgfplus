package msdbsearch.test;

import java.io.File;

import msdbsearch.SuffixArrayForMSGFDB;

import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;

public class SATest {
	public static void main(String argv[]) throws Exception
	{
		test();
	}
	
	public static void test() throws Exception
	{
		File databaseFile = new File("/Users/sangtaekim/Research/Data/IPI/IPI_human_3.79.fasta");
		SuffixArraySequence sequence = new SuffixArraySequence(databaseFile.getPath());
		
		// Normal suffix array read/write
		
		System.out.println("SuffixArray");
		long time = System.currentTimeMillis();
		SuffixArray sa = new SuffixArray(sequence);
		System.out.println("Time: " + (System.currentTimeMillis()-time));
		
		System.out.println("SuffixArrayForMSGFDB");
		time = System.currentTimeMillis();
		SuffixArrayForMSGFDB sa2 = new SuffixArrayForMSGFDB(sequence);
		System.out.println("Time: " + (System.currentTimeMillis()-time));
		System.out.println(sa.toString().equals(sa2.toString()));
		
		
	}

}
