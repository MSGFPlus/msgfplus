package edu.ucsd.msjava.msdbsearch.test;

import java.io.File;

import edu.ucsd.msjava.msdbsearch.SuffixArrayForMSGFDB;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.suffixarray.SuffixArray;
import edu.ucsd.msjava.suffixarray.SuffixArraySequence;


public class SATest {
	public static void main(String argv[]) throws Exception
	{
		test();
	}
	
	public static void test() throws Exception
	{
		File databaseFile = new File("/home/sangtaekim/Research/Data/IPI/IPI_human_3.79.fasta");
		SuffixArraySequence sequence = new SuffixArraySequence(databaseFile.getPath());
		
		// Normal suffix array read/write
		
		long time;
//		System.out.println("SuffixArray");
//		time = System.currentTimeMillis();
//		SuffixArray sa = new SuffixArray(sequence);
//		System.out.println("Time: " + (System.currentTimeMillis()-time));
		
		System.out.println("SuffixArrayForMSGFDB");
		time = System.currentTimeMillis();
		SuffixArrayForMSGFDB sa2 = new SuffixArrayForMSGFDB(sequence);
		System.out.println("Time: " + (System.currentTimeMillis()-time));
		int numCandidates = sa2.getNumCandidatePeptides(AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys(), (383.8754f-(float)Composition.H)*3-(float)Composition.H2O, new Tolerance(2.5f, false));
		System.out.println("NumCandidatePeptides: " + numCandidates);
		int length10 = sa2.getNumDistinctPeptides(10);
		System.out.println("NumUnique10: " + length10);
	} 

}
