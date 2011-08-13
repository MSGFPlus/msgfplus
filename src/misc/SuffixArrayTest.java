package misc;


import msdbsearch.CompactFastaSequence;
import msdbsearch.SuffixArrayForMSGFDB;
import suffixarray.SuffixArray;
import suffixarray.SuffixArraySequence;

public class SuffixArrayTest {
	public static void main(String argv[]) throws Exception
	{
		if(argv.length != 1)
		{
			System.out.println("java SuffixArrayTest *.fasta");
			System.exit(-1);
		}
		testSA(argv[0]);
	}
	
	public static void testSA(String fastaFile) throws Exception
	{
//	    String fastaFile = System.getProperty("user.home")+"/Research/Data/IPI/IPI_human_3.79.fasta";
//	    String fastaFile = System.getProperty("user.home")+"/Research/Data/IPI/IPI_human_3.79_shuffle.fasta";
//	    String fastaFile = System.getProperty("user.home")+"/Research/Data/ISBControl/Mix_7/ORBITRAP/Database/18mix.fasta";
//		String fastaFile = System.getProperty("user.home")+"/Research/Data/IPI/tiny.fasta";
		
	    long time = System.currentTimeMillis();
	    CompactFastaSequence sequence = new CompactFastaSequence(fastaFile);
	    System.out.println("-- Loading fasta file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
	    
//	    time = System.currentTimeMillis();
//	    SuffixArrayForMSGFDB sa = new SuffixArrayForMSGFDB(sequence);
//	    System.out.println("-- Loading SuffixArray file time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
//
//	    time = System.currentTimeMillis();
	 
//	    for(String m : sa.getAllMatchedStrings("abr"))
//	    	System.out.println(m);
//	    System.out.println("-- Searching time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
//	    
//	    System.out.println("SA: " + sa.toString());
//	    for(int l=1; l<50; l++)
//	    	System.out.println(l+"\t" + sa.getNumDistinctSeq(l));
//	    System.out.println(12+"\tNumDistinctPep\t" + sa.getNumDistinctSeq(12));
//	    float peptideMass = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys().getPeptide("").getMass();
//	    float peptideMass = 4567.5f;
//	    Tolerance tolerance = new Tolerance(1, true);
//	    for(int i=0; i<10; i++)
//	    {
//		    time = System.currentTimeMillis();
//		    System.out.println(sa.getNumCandidatePeptides(peptideMass, tolerance));
//		    System.out.println("-- Searching time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
//	    }
	    time = System.currentTimeMillis();
//	    sa.printAllPeptides(AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys(), 5, 50);
	    System.out.println("-- Searching time: " + (System.currentTimeMillis() - time)/1000.0 + "s");
	}
}
