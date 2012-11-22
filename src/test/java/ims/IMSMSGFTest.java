package ims;

import java.io.File;

import org.junit.Test;

import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGF;

public class IMSMSGFTest {
	@Test
	public void testIMSMSGF()
	{
		File baseDir = new File("/Users/kims336/Research/Data/IMS/Sarc_DTAs");
		File msgfInputFile = new File(baseDir.getPath()+File.separator+"test_msgfInput.txt");
		File specDir = baseDir;
		
		String[] argv = {"-i", msgfInputFile.getPath(), "-d", specDir.getPath(), 
				"-m", "1", "-inst", "2", "-e", "1", "-fixMod", "1", "-x", "0"};


		ParamManager paramManager = new ParamManager("MSGF", "Test", "Test", "java -Xmx2000M -cp MSGFDB.jar ui.MSGF");
		paramManager.addMSGFParams();

		if(argv.length == 0)
		{
			paramManager.printUsageInfo();
			return;
		}

		// Parse parameters
		paramManager.parseParams(argv);
		MSGF.runMSGF(paramManager);
	}
}
