package ims;

import java.io.File;
import java.net.URISyntaxException;

import org.junit.Test;

import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGF;

public class IMSMSGFTest {

	@Test
	public void testIMSMSGF() throws URISyntaxException {
		File msFile = new File(IMSMSGFTest.class.getClassLoader().getResource("test.mgf").toURI());
		File fastaFile = new File(IMSMSGFTest.class.getClassLoader().getResource("human-uniprot-contaminants.fasta").toURI());

		String[] argv = {"-db", fastaFile.getAbsolutePath(), "-i", msFile.getAbsolutePath(),
				"-m", "1", "-inst", "2", "-e", "1", "-fixMod", "1", "-x", "0"};
		ParamManager paramManager = new ParamManager("MSGF", "Test", "Test", "java -Xmx2000M -cp MSGFDB.jar ui.MSGF");
		paramManager.addMSGFParams();

		// Parse parameters
		paramManager.parseParams(argv);
		MSGF.runMSGF(paramManager);
	}
}
