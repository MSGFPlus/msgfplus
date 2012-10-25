package msgfplus;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;

public class TestMSGFPlus {
	@Test
	public void testAddFeatures()
	{
//		File specFile = new File(System.getProperty("user.home")+"/Research/Data/Viktor/test.ms2");
//		File dbFile = new File(System.getProperty("user.home")+"/Research/Data/Viktor/test.fasta");
//		
//		File tempOutputFile = null;
//		try {
//			tempOutputFile = File.createTempFile("__TestMSGFPlus", ".mzid");
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		tempOutputFile.deleteOnExit();
//		
//		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), "-addFeatures", "1", "-o", tempOutputFile.getPath()};
//		
//		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
//		paramManager.addMSGFPlusParams();
//		assertTrue(paramManager.parseParams(argv) == null);
//		
//		assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
	}
	
	@Test
	public void testMzIdGen()
	{
//		File specFile = new File("/Users/kims336/Research/Data/CPTAC/test.mgf");
//		File dbFile = new File("/Users/kims336/Research/Data/CPTAC/test.fasta");
//		File modFile = new File("/Users/kims336/Research/Data/CPTAC/Mods.txt");
//
//		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), "-mod", modFile.getPath()};
//		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
//		paramManager.addMSGFPlusParams();
//		assertTrue(paramManager.parseParams(argv) == null);
//		assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
	}
	
	@Test
	public void testUserScoringParams()
	{
//		File specFile = new File("/Users/kims336/Research/Data/CPTAC/test.mgf");
//		File dbFile = new File("/Users/kims336/Research/Data/CPTAC/test.fasta");
//		File modFile = new File("/Users/kims336/Research/Data/CPTAC/Mods.txt");
//
//		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), "-mod", modFile.getPath(), "-protocol", "2", "-m", "3"};
//		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
//		paramManager.addMSGFPlusParams();
//		assertTrue(paramManager.parseParams(argv) == null);
//		assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
	}
}
