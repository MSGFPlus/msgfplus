package ims;

import java.io.File;

import org.junit.Test;

import edu.ucsd.msjava.msdbsearch.PeptideEnumerator;

public class SarcTest {
	@Test
	public void generatePeptideList()
	{
		File fastaFile = new File("C:\\cygwin\\home\\kims336\\Data\\IMS_Sarc\\IPI_human_3.87_withContam_reverse.fasta");
//		File fastaFile = new File("C:\\cygwin\\home\\kims336\\Data\\IMS_BSA\\BSA.fasta");
//		File outputFile = new File("C:\\cygwin\\home\\kims336\\Data\\IMS_BSA\\BSAPeptides_ST.txt");
		File outputFile = new File("C:\\cygwin\\home\\kims336\\Data\\IMS_Sarc\\HumanPeptides_Reverse.txt");
		try {
			PeptideEnumerator.enumerate(fastaFile, outputFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("Done");
	}
}
