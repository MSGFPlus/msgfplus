package msgfplus;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import edu.ucsd.msjava.mzml.MzMLAdapter;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;
import edu.ucsd.msjava.ui.ScoringParamGen;

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
	
	@Test
	public void testMultipleInputSpectra()
	{
//		File specPath = new File("/Users/kims336/Research/Data/TrainingITRAQ/Global/mzMLFiles");
//		File dbFile = new File("/Users/kims336/Research/Data/TrainingITRAQ/Global/ID_003632_9011437E.fasta");
//		File modFile = new File("/Users/kims336/Research/Data/TrainingITRAQ/Global/Mods.txt");
//
//		String[] argv = {"-s", specPath.getPath(), "-d", dbFile.getPath(), "-mod", modFile.getPath(), "-ti", "0,1", "-tda", "1", "-ntt", "2", "-protocol", "2"};
//		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
//		paramManager.addMSGFPlusParams();
//		paramManager.parseParams(argv);
//		MSGFPlus.runMSGFPlus(paramManager);
	}
	
	@Test
	public void testNTermMetCleavage()
	{
		File specPath = new File("C:\\cygwin\\home\\kims336\\Research\\Data\\Bug\\test.mgf");
		File dbFile = new File("C:\\cygwin\\home\\kims336\\Research\\Data\\Bug\\test.fasta");

		String[] argv = {"-s", specPath.getPath(), "-d", dbFile.getPath(), "-ti", "0,0", "-tda", "0", "-ntt", "2", "-protocol", "2"};
		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		paramManager.parseParams(argv);
		MSGFPlus.runMSGFPlus(paramManager);
	}
	
	@Test
	public void testScoringParamGen()
	{
		File resultPath = new File("/Users/kims336/Research/Data/TrainingITRAQ/Phospho/test");
		File specPath = new File("/Users/kims336/Research/Data/TrainingITRAQ/Phospho/mzMLFiles");
		
		String[] argv = {"-i", resultPath.getPath(), "-d", specPath.getPath(), "-m", "2", "-inst", "1", "-e", "0", "-protocol", "3"};
		ParamManager paramManager = new ParamManager("ScoringParamGen", String.valueOf(ScoringParamGen.VERSION), ScoringParamGen.DATE,
				"java -Xmx2000M -cp MSGFPlus.jar edu.ucsd.msjava.ui.ScoringParamGen");
		MzMLAdapter.turnOffLogs();
		paramManager.addScoringParamGenParams();
		paramManager.parseParams(argv);
		
		ScoringParamGen.runScoringParamGen(paramManager);
	}
	
	@Test
	public void testMgf()
	{
		File specPath = new File("C:\\cygwin\\home\\kims336\\Research\\Data\\Test\\test.mgf");
		File dbFile = new File("C:\\cygwin\\home\\kims336\\Research\\Data\\Test\\test.fasta");

		String[] argv = {"-s", specPath.getPath(), "-d", dbFile.getPath(), "-ti", "0,0", "-tda", "0", "-ntt", "2"};
		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		paramManager.parseParams(argv);
		MSGFPlus.runMSGFPlus(paramManager);
	}

}
