package msgfplus;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import edu.ucsd.msjava.mzid.MzIDParser;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;

public class TestMSGFPlus {
	@Test
	public void testAddFeatures()
	{
		File specFile = new File(System.getProperty("user.home")+"/Research/Data/Viktor/test.ms2");
		File dbFile = new File(System.getProperty("user.home")+"/Research/Data/Viktor/test.fasta");
		
		File tempOutputFile = null;
		try {
			tempOutputFile = File.createTempFile("__TestMSGFPlus", ".mzid");
		} catch (IOException e) {
			e.printStackTrace();
		}
		tempOutputFile.deleteOnExit();
		
		System.out.println(tempOutputFile.getPath());
		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), "-addFeatures", "1"};//, "-o", tempOutputFile.getPath()};
		
		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		assertTrue(paramManager.parseParams(argv) == null);
		
		assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
		
//		MzIDParser parser = new MzIDParser(tempOutputFile);
	}
}
