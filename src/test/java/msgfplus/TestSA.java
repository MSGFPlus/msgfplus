package msgfplus;

import java.io.File;

import org.junit.Test;

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

}
