package msgfplus;

import static org.junit.Assert.*;

import java.io.File;
import java.util.HashSet;
import java.util.Iterator;

import org.junit.Test;

import edu.ucsd.msjava.msgf.Histogram;
import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.IonType;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.mzid.MzIDTest;
import edu.ucsd.msjava.mzml.MzMLAdapter;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.ui.MSGFPlus;
import edu.ucsd.msjava.ui.ScoringParamGen;

public class TestMSGFPlus {
	
	@Test
	public void testSingleSpec()
	{
		File dir = new File("C:\\cygwin\\home\\kims336\\Data\\Debug");

		File specFile = new File(dir.getPath()+File.separator+"test.mgf");
		File dbFile = new File(dir.getPath()+File.separator+"test.fasta");
		File modFile = new File(dir.getPath()+File.separator+"Mods.txt");
		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), 
				"-mod", modFile.getPath(), "-t", "25ppm", "-ti", "0,1"//, "-tda", "1", "-m", "1"
//				, "-maxLength", "250"
				}; 

		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		
		String msg = paramManager.parseParams(argv);
		assertTrue(msg == null);
		
		assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
	}

	@Test
	public void testQCShew()
	{
		File dir = new File("C:\\cygwin\\home\\kims336\\Data\\QCShew");

		File specFile = new File(dir.getPath()+File.separator+"QC_Shew_12_02_2_1Aug12_Cougar_12-06-11_dta.txt");
		File dbFile = new File(dir.getPath()+File.separator+"ID_003456_9B916A8B.fasta");
		File modFile = new File(dir.getPath()+File.separator+"Mods.txt");
//		File outputFile = new File(dir.getPath()+File.separator+"Test"+"2013-07-26"+".txt");
		String versionString = MSGFPlus.VERSION.split("\\s+")[1];
		versionString = versionString.substring(versionString.indexOf('(')+1, versionString.lastIndexOf(')'));
		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), 
				"-mod", modFile.getPath(), "-t", "10ppm", "-tda", "1", "-m", "1", "-ti", "0,1", "-ntt", "1",
				"-o", dir.getPath()+File.separator+"Test_"+versionString+".mzid"
				}; 

		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		
		String msg = paramManager.parseParams(argv);
		if(msg != null)
			System.err.println("Error: " + msg);
		assertTrue(msg == null);
		
		assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
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
//		File specFile = new File("C:\\cygwin\\home\\kims336\\Data\\Centroiding\\QC_Shew_13_02_pt5_2_pH7p5_23Apr13_Frodo_12-12-04.mgf");
		File specFile = new File("C:\\cygwin\\home\\kims336\\Data\\Centroiding\\H20120525_JQ_CPTAC2_Compref4_protfxn01.mgf");
		SpectraAccessor specAcc = new SpectraAccessor(specFile);
		Iterator<Spectrum> itr = specAcc.getSpecItr();
		int numSpecs = 0;
		int numProfileSpecs = 0;
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			numSpecs++;
			spec.determineIsCentroided();
			if(!spec.isCentroided())
				numProfileSpecs++;
		}
		System.out.println("NumSpecs: " + numSpecs);
		System.out.println("NumProfileSpecs: " + numProfileSpecs);
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
	
	@Test
	public void testPrePost() throws Exception
	{
		File specPath = new File("C:\\cygwin\\home\\kims336\\Research\\Data\\Matt\\test.mgf");
		File dbFile = new File("C:\\cygwin\\home\\kims336\\Research\\Data\\Matt\\test.fasta");
		File modFile = new File("C:\\cygwin\\home\\kims336\\Research\\Data\\Matt\\MSGFDB_Mods.txt");

		String[] argv = {"-s", specPath.getPath(), "-d", dbFile.getPath(), "-t", "1.5Da,2.5Da", "-n", "2", "-ntt", "1", "-mod", modFile.getPath(), "-ti", "0,1"};
		
		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		paramManager.parseParams(argv);
		MSGFPlus.runMSGFPlus(paramManager);
		System.out.println("Done");
	}

	@Test
	public void testNoEnzymeSearch() throws Exception
	{
		File specPath = new File("C:\\cygwin\\home\\kims336\\Research\\Data\\Matt\\test.mgf");
		File dbFile = new File("C:\\cygwin\\home\\kims336\\Research\\Data\\Matt\\test.fasta");
		File modFile = new File("C:\\cygwin\\home\\kims336\\Research\\Data\\Matt\\MSGFDB_Mods.txt");

		String[] argv = {"-s", specPath.getPath(), "-d", dbFile.getPath(), "-t", "1.5Da,2.5Da", "-n", "2", "-ntt", "1", "-mod", modFile.getPath(), "-ti", "0,1", "-e", "0"};
		
		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		paramManager.parseParams(argv);
		MSGFPlus.runMSGFPlus(paramManager);
		System.out.println("Done");
	}

	@Test
	public void testCTermFixedMod() throws Exception
	{
		File specPath = new File(System.getProperty("user.home")+"/Research/Data/Joe/test.mgf");
		File dbFile = new File(System.getProperty("user.home")+"/Research/Data/Joe/test.fasta");
		File modFile = new File(System.getProperty("user.home")+"/Research/Data/Joe/testMods.txt");

//		String[] argv = {"-s", specPath.getPath(), "-d", dbFile.getPath(), "-t", "2.5Da,2.5Da", "-n", "1", "-ntt", "2", "-mod", modFile.getPath(), "-ti", "1"};
		String[] argv = {"-s", specPath.getPath(), "-d", dbFile.getPath(), "-t", "10ppm", "-n", "1", "-ntt", "2", "-mod", modFile.getPath(), "-ti", "1"};
		
		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		paramManager.parseParams(argv);
		MSGFPlus.runMSGFPlus(paramManager);
		System.out.println("Done");
	}
	
	@Test
	public void testScoringParamGenETD()
	{
		File resultPath = new File("/Users/kims336/Research/Data/TrainingITRAQ/Global/beforeTraining/CPTAC_OvC_JB5427_iTRAQ_01_9Apr12_Cougar_12-03-21.mzid");
		File specPath = new File("/Users/kims336/Research/Data/TrainingITRAQ/Global/mzMLFiles");
		
		String[] argv = {"-i", resultPath.getPath(), "-d", specPath.getPath(), "-m", "1", "-inst", "0", "-e", "1"};
		ParamManager paramManager = new ParamManager("ScoringParamGen", String.valueOf(ScoringParamGen.VERSION), ScoringParamGen.DATE,
				"java -Xmx2000M -cp MSGFPlus.jar edu.ucsd.msjava.ui.ScoringParamGen");
		MzMLAdapter.turnOffLogs();
		paramManager.addScoringParamGenParams();
		paramManager.parseParams(argv);
		
		ScoringParamGen.runScoringParamGen(paramManager);
	}

	@Test
	public void selectIonTypeTest()
	{
		IonType.getAllKnownIonTypes(4, true);
	}
	
	@Test
	public void testQExactiveParam()
	{
		File specPath = new File(System.getProperty("user.home")+"/Research/Data/Tabb/TCGA-7D-HCD-QExactive/test.mgf");
		File dbFile = new File(System.getProperty("user.home")+"/Research/Data/Tabb/TCGA-7D-HCD-QExactive/test.fasta");

		String[] argv = {"-s", specPath.getPath(), "-d", dbFile.getPath(), "-t", "10ppm", "-ntt", "1", "-ti", "0,0", "-m", "3", "-inst", "3"};
		
		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		paramManager.parseParams(argv);
		MSGFPlus.runMSGFPlus(paramManager);
		System.out.println("Done");		
	}
	
	@Test
	public void testScoringProtocol()
	{
		File specPath = new File(System.getProperty("user.home")+"/Research/Data/Tabb/TCGA-7D-HCD-QExactive/test.mgf");
		File dbFile = new File(System.getProperty("user.home")+"/Research/Data/Tabb/TCGA-7D-HCD-QExactive/test.fasta");
		File modFile = new File(System.getProperty("user.home")+"/Research/Data/Tabb/TCGA-7D-HCD-QExactive/testMod.txt");

		String[] argv = {"-s", specPath.getPath(), "-d", dbFile.getPath(), "-t", "10ppm", "-ntt", "1", "-ti", "0,0", "-m", "3", "-inst", "3", "-mod", modFile.getPath()};
		
		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		paramManager.parseParams(argv);
		MSGFPlus.runMSGFPlus(paramManager);
		System.out.println("Done");		
	}	
	
	@Test
	public void testMultipleNTermMod()
	{
		File specPath = new File(System.getProperty("user.home")+"/Research/Data/Debug/test.mgf");
		File dbFile = new File(System.getProperty("user.home")+"/Research/Data/Debug/test.fasta");
		File modFile = new File(System.getProperty("user.home")+"/Research/Data/Debug/mods.txt");

		String[] argv = {"-s", specPath.getPath(), "-d", dbFile.getPath(), "-t", "20ppm", "-ntt", "2", "-ti", "0,1", "-m", "3", "-inst", "1", "-mod", modFile.getPath()};
		
		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		paramManager.parseParams(argv);
		MSGFPlus.runMSGFPlus(paramManager);
		System.out.println("Done");		
	}
	
	@Test
	public void testTMT()
	{
		File specPath = new File(System.getProperty("user.home")+"/Research/Data/Tabb/TCGA-7D-HCD-QExactive/test.mgf");
		File dbFile = new File(System.getProperty("user.home")+"/Research/Data/Tabb/TCGA-7D-HCD-QExactive/test.fasta");
		File modFile = new File(System.getProperty("user.home")+"/Research/Data/Tabb/TCGA-7D-HCD-QExactive/testMod.txt");

		String[] argv = {"-s", specPath.getPath(), "-d", dbFile.getPath(), "-t", "10ppm", "-ntt", "1", "-ti", "0,0", "-m", "3", "-inst", "3", "-mod", modFile.getPath()};
		
		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		paramManager.parseParams(argv);
		MSGFPlus.runMSGFPlus(paramManager);
		System.out.println("Done");		
	}		
	
	@Test
	public void testWaLP()
	{
//		File dir = new File("/Users/kims336/Research/Data/Jesse");
//		File specPath = new File(dir.getPath() + File.separator + "test.mgf");
//		File dbFile = new File(dir.getPath() + File.separator + "test.fasta");

		File specPath = new File(System.getProperty("user.home")+"/Research/Data/Jesse/test.mgf");
		File dbFile = new File(System.getProperty("user.home")+"/Research/Data/Jesse/test.fasta");
		
		String[] argv = {"-s", specPath.getPath(), "-d", dbFile.getPath(), "-t", "0.1Da", "-m", "1", "-inst", "1", "-e", "10"};
		
		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		paramManager.parseParams(argv);
		MSGFPlus.runMSGFPlus(paramManager);
		System.out.println("Done");		
	}	
	
	@Test
	public void testSharedPeptides()
	{
		File dir = new File("D:\\Research\\Data\\MSGFPercolator");

		File specFile = new File(dir.getPath()+File.separator+"testCID.mgf");
		File dbFile = new File(dir.getPath()+File.separator+"testCID.fasta");
		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath()};

//		File specFile = new File(dir.getPath()+File.separator+"testITRAQ.mgf");
//		File dbFile = new File(dir.getPath()+File.separator+"testITRAQ.fasta");
//		File modFile = new File(dir.getPath()+File.separator+"testITRAQMods.txt");
//		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), "-addFeatures", "1", "-mod", modFile.getPath(), "-m", "3"};
		
		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		
		String msg = paramManager.parseParams(argv);
		assertTrue(msg == null);
		
		assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
	}	
	
	@Test
	public void testAdditionalFeatures()
	{
		File dir = new File("C:\\cygwin\\home\\sangtaekim\\Research\\Data\\CompRef\\Global\\");

		File specFile = new File(dir.getPath()+File.separator+"test.mgf");
		File dbFile = new File(dir.getPath()+File.separator+"test.fasta");
		File modFile = new File(dir.getPath()+File.separator+"Mods.txt");
		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), 
				"-mod", modFile.getPath(), "-t", "10ppm", "-m", "3", //}; 
				"-addFeatures", "1"};

		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		
		String msg = paramManager.parseParams(argv);
		assertTrue(msg == null);
		
		assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
		
	}
	
	@Test
	public void testLongPeptides()
	{
		File dir = new File("C:\\cygwin\\home\\kims336\\Data\\QCShew");

		File specFile = new File(dir.getPath()+File.separator+"QC_Shew_12_02_pt5_a4_10Apr13_Draco_13-02-14_dta.txt");
		File dbFile = new File(dir.getPath()+File.separator+"ID_003456_9B916A8B.fasta");
		File modFile = new File(dir.getPath()+File.separator+"Mods.txt");
		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), 
				"-mod", modFile.getPath(), "-t", "10ppm", "-tda", "1", "-m", "3",
				"-maxLength", "250"}; 

		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		
		String msg = paramManager.parseParams(argv);
		assertTrue(msg == null);
		
		assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
	}
	
	@Test
	public void testCentroiding()
	{
		File specFile = new File("C:\\cygwin\\home\\kims336\\Data\\Debug\\ALZ_VP2P101_C_SCX_02_7Dec08_Draco_08-10-29_dta.txt");
		SpectraAccessor acc = new SpectraAccessor(specFile);
		while(acc.getSpecItr().hasNext())
		{
			Spectrum spec = acc.getSpecItr().next();
			spec.determineIsCentroided();
			if(!spec.isCentroided())
			{
				System.out.println("ScanNum: " + spec.getScanNum());
				System.exit(-1);
			}
		}
	}
	
	@Test
	public void testNumVariantsPerPeptide()
	{
		File dir = new File("C:\\cygwin\\home\\kims336\\Data\\Debug");

		File specFile = new File(dir.getPath()+File.separator+"test.mgf");
		File dbFile = new File(dir.getPath()+File.separator+"test.fasta");
		File modFile = new File(dir.getPath()+File.separator+"Mods.txt");
		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), 
				"-mod", modFile.getPath(), "-t", "10ppm", "-tda", "1", "-m", "1", "-ti", "0,1"
//					, "-maxLength", "250"
				}; 

		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		
		String msg = paramManager.parseParams(argv);
		assertTrue(msg == null);
		
		assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
	}
	
	@Test
	public void charge1Hisgogram()
	{
		File specFile = new File("D:\\Research\\Data\\QCShew\\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11_mzML.mgf");
		SpectraAccessor specAcc = new SpectraAccessor(specFile);
		Histogram<Integer> ticHist = new Histogram<Integer>();
		while(specAcc.getSpecItr().hasNext())
		{
			Spectrum spec = specAcc.getSpecItr().next();
			int charge = spec.getCharge();

			if(charge == 1)
			{
				float tic = 0;
				float ticBelowPrecursor = 0;
				float precursorMz = spec.getPrecursorPeak().getMz();
				for(Peak p : spec)
				{
					tic += p.getIntensity();
					if(p.getMz() < precursorMz)
						ticBelowPrecursor += p.getIntensity();
				}
				ticHist.add(Math.round(ticBelowPrecursor/tic*100));
			}
		}
		ticHist.printSortedRatio();
	}
	
	@Test
	public void testCharge1Detection()
	{
		File dir = new File("D:\\Research\\Data\\QCShew\\");

		File specFile = new File(dir.getPath()+File.separator+"charge1.mgf");
		File dbFile = new File(dir.getPath()+File.separator+"charge1.fasta");
		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), 
				"-t", "20ppm"
				}; 

		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		
		String msg = paramManager.parseParams(argv);
		assertTrue(msg == null);
		
		assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);		
	}
	
	@Test
	public void testOrbiHCDScoring()
	{
		File dir = new File("C:\\cygwin\\home\\kims336\\Data\\Debug");

		File specFile = new File(dir.getPath()+File.separator+"test.mgf");
		File dbFile = new File(dir.getPath()+File.separator+"test.fasta");
		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), 
				"-m", "3", "-ntt", "1", "-inst", "1", "-n", "1", "-addFeatures", "1"
				}; 

		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		
		String msg = paramManager.parseParams(argv);
		if(msg != null)
			System.out.println(msg);
		assertTrue(msg == null);
		
		assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
	}

	@Test
	public void testCCMSXMLParsing()
	{
		File dir = new File("C:\\cygwin\\home\\kims336\\Data\\Debug");
		File paramFile = new File(dir.getPath()+File.separator+"params.xml");
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromXMLFile(paramFile.getPath());
		aaSet.printAASet();
	}
	
	@Test
	public void testDoubleMods()
	{
		File mod = new File("C:\\cygwin\\home\\kims336\\Data\\Debug\\Mods.txt");
		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile(mod.getPath());
		aaSet.printAASet();
	}
	
	@Test
	public void testMzIdOutput()
	{
//		File mzidFile = new File("C:\\cygwin\\home\\kims336\\Data\\Mayo\\reproducibilitySample.mzid");
		File mzidFile = new File("C:\\cygwin\\home\\kims336\\Data\\Mayo\\AllAll.mzid");
		MzIDTest test = new MzIDTest(mzidFile);
		boolean isValid;
		try {
			isValid = test.isValid();
			System.out.println("IsValid? " + isValid);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	@Test
	public void testDuplicatePepEv()
	{
		File dir = new File("C:\\cygwin\\home\\kims336\\Data\\Mayo");

//		File specFile = new File(dir.getPath()+File.separator+"test.mgf");
		File specFile = new File(dir.getPath()+File.separator+"reproducibilitySample.mgf");
		File dbFile = new File(dir.getPath()+File.separator+"H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta");
		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), 
				"-n", "5", "-o", dir.getPath()+File.separator+"AllAll.mzid"
				}; 

		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		
		String msg = paramManager.parseParams(argv);
		assertTrue(msg == null);
		
		assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
	}
	
}
