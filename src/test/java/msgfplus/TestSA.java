package msgfplus;

import java.io.File;

import org.junit.Test;

import edu.ucsd.msjava.msdbsearch.CompactFastaSequence;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.suffixarray.SuffixArray;
import edu.ucsd.msjava.suffixarray.SuffixArraySequence;

public class TestSA {
	@Test
	public void getNumCandidatePeptides()
	{
		File dbFile = new File("/Users/kims336/Research/Data/Andy/IPI_human_3.87_withContam.fasta");
		SuffixArraySequence sequence = new SuffixArraySequence(dbFile.getPath());
		SuffixArray sa = new SuffixArray(sequence);
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile(System.getProperty("user.home")+"/Research/Data/Andy/TestMods.txt");
		System.out.println("NumPeptides: " + sa.getNumCandidatePeptides(aaSet, 2364.981689453125f, new Tolerance(10, true)));
	}

	@Test
	public void testRedundantProteins()
	{
		File databaseFile = new File("D:\\Research\\Data\\IPRG2014\\sprot-ecoli-4spiked-20131016.fasta");
		
		CompactFastaSequence fastaSequence = new CompactFastaSequence(databaseFile.getPath());
		float ratioUniqueProteins = fastaSequence.getRatioUniqueProteins();
		if(ratioUniqueProteins < 0.5f)
		{
			System.err.println("Error while indexing: " + databaseFile.getName() + " (too many redundant proteins)");
			System.err.println("UniqueProteinRation: " + ratioUniqueProteins);
			System.err.println("If the database contains forward and reverse proteins, run MS-GF+ (or BuildSA) again with \"-tda 0\"");
			System.exit(-1);
		}
		
		float fractionDecoyProteins = fastaSequence.getFractionDecoyProteins();
		if(fractionDecoyProteins < 0.4f || fractionDecoyProteins > 0.6f)
		{
			System.err.println("Error while reading: " + databaseFile.getName() + " (fraction of decoy proteins: "+ fractionDecoyProteins+ ")");
			System.err.println("Delete " + databaseFile.getName() + " and run MS-GF+ again.");
			System.exit(-1);
		}
		
	}
	
}
