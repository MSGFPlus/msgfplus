package ims;

import java.io.File;

import edu.ucsd.msjava.ui.MSGFPlus;
import org.junit.Test;
import edu.ucsd.msjava.msdbsearch.PeptideEnumerator;

public class SarcTest {
	@Test
	public void generatePeptideList() throws Exception {
		File fastaFile = new File(SarcTest.class.getClassLoader().getResource("human-uniprot-contaminants.fasta").toURI());
		File outputFile = File.createTempFile("HumanPeptides", "tsv");
		PeptideEnumerator.enumerate(fastaFile, outputFile, MSGFPlus.DEFAULT_DECOY_PROTEIN_PREFIX);
		outputFile.deleteOnExit();
		System.out.println("Done");
	}
}
