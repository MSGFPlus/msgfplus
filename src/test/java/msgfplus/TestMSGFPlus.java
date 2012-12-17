package msgfplus;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;

import org.junit.Test;

import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.mzml.MzMLAdapter;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.parser.BufferedLineReader;
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

	@Test
	public void testCrossLink()
	{
		File specPath = new File(System.getProperty("user.home")+"/Research/Data/Eric/CuMB_006e_scan2089.mgf");
		File dbFile = new File(System.getProperty("user.home")+"/Research/Data/Eric/database.fasta");

		String[] argv = {"-s", specPath.getPath(), "-d", dbFile.getPath(), "-ti", "0,0", "-tda", "0", "-ntt", "2", "-m", "2"};
		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		paramManager.parseParams(argv);
		MSGFPlus.runMSGFPlus(paramManager);
	}
	
	@Test
	public void testModFileReading()
	{
		File modFile = new File("/Users/kims336/Research/Data/Mouse_Brain_Phospho/Mods.txt");
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile.getPath());
		aaSet.printAASet();
	}
	
	@Test
	public void testProfieSpectraDetector()
	{
		File specFile = new File("/Users/kims336/Research/Data/Nikola/Athal0470_26Mar12_Jaguar_12-02-27_dta.txt");
		SpectraAccessor specAcc = new SpectraAccessor(specFile);
		Spectrum spec = specAcc.getSpectrumBySpecIndex(1);
		spec.setIsCentroided();
		System.out.println(spec.isCentroided());
	}
	
	@Test
	public void testDuplicatedProteins() throws Exception
	{
		File tsvFile = new File("/Users/kims336/Research/Data/Nikola/NoEnzyme.tsv");
		BufferedLineReader in = new BufferedLineReader(tsvFile.getPath());
		
		String headerLine = in.readLine();
		String[] header = headerLine.split("\t");
		int proteinCol = -1;
		for(int i=0; i<header.length; i++)
			if(header[i].equals("Protein"))
				proteinCol = i;
		
		String s;
		while((s=in.readLine()) != null)
		{
			if(s.startsWith("#"))
				continue;
			String[] token = s.split("\t");
			String annotation = token[proteinCol];
			String[] protein = annotation.split(";");
			HashSet<String> proteinSet = new HashSet<String>();
			for(String prot : protein)
			{
				if(!proteinSet.add(prot))
				{
					System.err.println("Duplicate Protein entry (" + prot + ")\n");
					System.err.println(s);
					System.exit(-1);
				}
			}
		}
		
		in.close();
		
		System.out.println("No duplicate protein entry");
	}
	
	@Test
	public void testPrecursorError() throws Exception
	{
		File specPath = new File(System.getProperty("user.home")+"/Research/Data/Matt/test_dta.txt");
		File dbFile = new File(System.getProperty("user.home")+"/Research/Data/Matt/test.fasta");
		File modFile = new File(System.getProperty("user.home")+"/Research/Data/Matt/Mods.txt");

		String[] argv = {"-s", specPath.getPath(), "-d", dbFile.getPath(), "-t", "1.5Da,2.5Da", "-ti", "-1,1", "-tda", "0", "-ntt", "1", "-n", "3", "-mod", modFile.getPath()};
		
		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		paramManager.parseParams(argv);
		MSGFPlus.runMSGFPlus(paramManager);
		System.out.println("Done");
	}
	
	@Test
	public void testMissingPeptides() throws Exception
	{
		File specPath = new File(System.getProperty("user.home")+"/Research/Data/Matt/test.mgf");
//		File dbFile = new File(System.getProperty("user.home")+"/Research/Data/Matt/ID_003456_9B916A8B.fasta");
		File dbFile = new File(System.getProperty("user.home")+"/Research/Data/Matt/test.fasta");

		String[] argv = {"-s", specPath.getPath(), "-d", dbFile.getPath(), "-t", "1.5Da,2.5Da", "-n", "1", "-ntt", "2"};
		
		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		paramManager.parseParams(argv);
		MSGFPlus.runMSGFPlus(paramManager);
		System.out.println("Done");
	}
	
}
