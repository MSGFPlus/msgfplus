package msgfplus;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

import edu.ucsd.msjava.ipa.IPA;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.parser.BufferedLineReader;
import edu.ucsd.msjava.parser.TSVResultParser;
import edu.ucsd.msjava.ui.MSGFPlus;

public class TestIPA {
	@Test
	public void testIPA() throws Exception
	{
		File dir = new File("C:\\cygwin\\home\\sangtaekim\\Research\\Data\\ASMS2013\\IPA\\");
		File peaksFile = new File(dir.getPath() + "\\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11_peaks.txt");
		File msgfPlusTsvFile = new File(dir.getPath() + "\\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11_NoCharge.tsv");

		File ipaResultFile = new File(dir.getPath() + "\\IPA.tsv");
		
		IPA ipa = new IPA(peaksFile, msgfPlusTsvFile);
		ipa.writeTo(ipaResultFile);
		System.out.println("Done");
	}
	
	@Test
	public void testReadingPeaks() throws Exception
	{
		File peaksFile = new File("D:\\Research\\Data\\ASMS2013\\NewDeconMSn\\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11_peaks.txt");
		BufferedLineReader in = new BufferedLineReader(peaksFile.getPath());
		String s;
		int lineNum = 0;
		String[] token;
		in.readLine();	// header
		long time = System.currentTimeMillis();
		while((s=in.readLine()) != null)
		{
			++lineNum;
			token = s.split("\t");
			int scanNum = Integer.parseInt(token[1]);
			double mz = Double.parseDouble(token[2]);
			double intensity = Double.parseDouble(token[3]);
		}
		System.out.println(lineNum);
		System.out.println("Elapsed Time: " + (System.currentTimeMillis() - time));
	}

	@Test
	public void compareIPAAndMSGFPlus()
	{
		File dir = new File("C:\\cygwin\\home\\sangtaekim\\Research\\Data\\ASMS2013\\");

		File ipaResultFile = new File(dir.getPath() + "\\IPA\\IPA_OnePerScan.tsv");
		TSVResultParser ipaParser = new TSVResultParser(ipaResultFile);
		ipaParser.parse(0.01f);
		Set<String> pepSet = ipaParser.getPepSet();
		System.out.println("IPAPeptides: " + pepSet.size());

		File msgfPlusResultFile = new File(dir.getPath() + "\\OldDeconMSn\\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11_dta_OnePerScan.tsv");
		TSVResultParser msgfplusParser = new TSVResultParser(msgfPlusResultFile);
		msgfplusParser.parse(0.01f);
		
		int numMSGFPlusOnly = 0;
		for(String pep : msgfplusParser.getPepSet())
		{
			if(!pepSet.contains(pep))
			{
				++numMSGFPlusOnly;
				System.out.println(pep);
			}
		}
		System.out.println("MSGFPlusOnly: " + numMSGFPlusOnly);
	}
	
	@Test
	public void compareMaxQuantAndMSGFPlus()
	{
		File dir = new File("C:\\cygwin\\home\\sangtaekim\\Research\\Data\\ASMS2013\\");
//		File msgfPlusResultFile = new File(dir.getPath() + "\\OldDeconMSn\\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11_dta.tsv");
		File msgfPlusResultFile = new File(dir.getPath() + "\\IPA\\IPA_OnePerScan.tsv");
		TSVResultParser msConvertResParser = new TSVResultParser(msgfPlusResultFile);
		msConvertResParser.parse(0.01f);
		Set<String> pepSet = msConvertResParser.getPepSet();
//		for(String pep : pepSet)
//		{
//			System.out.println(pep);
//		}
		System.out.println("NumMSGFPlusPeptides: " + pepSet.size());

		int numMQOnly = 0;
		File maxQuantResultFile = new File(dir.getPath() + "\\Shewanella\\combined\\txt\\peptides.txt");
		try {
			BufferedLineReader in = new BufferedLineReader(maxQuantResultFile.getPath());
			String s;
			
			in.readLine();	// header
			while((s=in.readLine()) != null)
			{
				String[] token = s.split("\\s+");
				String mqPep = token[0];
				if(!pepSet.contains(mqPep))
				{
					numMQOnly++;
					System.out.println(mqPep);
				}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		System.out.println("MaxQuantOnly: " + numMQOnly);
	}
	
	@Test
	public void iPRG2013()
	{
		File dir = new File("D:\\Research\\Data\\iPRG2013");

		File specFile = new File(dir.getPath()+File.separator+"DeconMSn\\f01.mgf");
		File dbFile = new File(dir.getPath()+File.separator+"Homo_sapiens_non-redundant.GRCh37.68.pep.all_FPKM_NOVEL-cRAP_targetdecoy.fasta");
//		File specFile = new File(dir.getPath()+File.separator+"test.mgf");
//		File dbFile = new File(dir.getPath()+File.separator+"test.fasta");
		File modFile = new File(dir.getPath()+File.separator+"Mods.txt");

		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), 
				"-mod", modFile.getPath(), "-t", "10ppm", "-tda", "0", "-ti", "0,0", "-inst" , "1", "-protocol", "4", "-m", "3"
		}; 

		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		
		String msg = paramManager.parseParams(argv);
		assertTrue(msg == null);
		
		assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
		
	}

	@Test
	public void countDiffScans()
	{
		File dir = new File("D:\\Research\\Data\\iPRG2013");
		File deconMSnSpecFile = new File(dir.getPath() + File.separator + "DeconMSn0515\\f01.mgf");
		SpectraAccessor deconMSnSpecAccessor = new SpectraAccessor(deconMSnSpecFile);
		Iterator<Spectrum> itr = deconMSnSpecAccessor.getSpecItr();
		
		int numSpectra = 0;
		HashSet<Integer> scanSet = new HashSet<Integer>();
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			++numSpectra;
			if(scanSet.add(spec.getScanNum()) == true)
			{
				System.out.println(spec.getScanNum());
			}
		}
		System.out.println("NumSpectra: " + numSpectra);
		System.out.println("NumScans: " + scanSet.size());
	}
	
	@Test
	public void compareIPRG2013Results()
	{
		File dir = new File("D:\\Research\\Data\\iPRG2013");

		File msConvertSpecFile = new File(dir.getPath() + File.separator + "mzML\\F01.mzML");
		SpectraAccessor v1SpecAccessor = new SpectraAccessor(msConvertSpecFile);
		
		File msConvertResultFile = new File(dir.getPath() + File.separator + "PW_F01_TI0_FDR.tsv");
		TSVResultParser msConvertResParser = new TSVResultParser(msConvertResultFile);
		msConvertResParser.parse(0.01f);
		Set<String> msConvertScans = msConvertResParser.getScanSet();
		Map<String,Integer> idToScan1 = new HashMap<String,Integer>();
		Map<Integer,String> scanToId1 = new HashMap<Integer,String>();
		Iterator<Spectrum> v1SpecItr = v1SpecAccessor.getSpecItr();
		while(v1SpecItr.hasNext())
		{
			Spectrum spec = v1SpecItr.next();
			idToScan1.put(spec.getID(), spec.getScanNum());
			scanToId1.put(spec.getScanNum(), spec.getID());
		}
		
		File deconMSnSpecFile = new File(dir.getPath() + File.separator + "DeconMSn0516\\F01.mgf");
		SpectraAccessor deconMSnSpecAccessor = new SpectraAccessor(deconMSnSpecFile);

		Map<String,Integer> idToScan2 = new HashMap<String,Integer>();
		Map<Integer,String> scanToId2 = new HashMap<Integer,String>();
		Iterator<Spectrum> v2SpecItr = deconMSnSpecAccessor.getSpecItr();
		while(v2SpecItr.hasNext())
		{
			Spectrum spec = v2SpecItr.next();
			idToScan2.put(spec.getID(), spec.getScanNum());
			scanToId2.put(spec.getScanNum(), spec.getID());
		}
		
		File deconMSnResultFile = new File(dir.getPath() + File.separator + "DeconMSn0516\\F01_FDR.tsv");
		TSVResultParser deconMSnResParser = new TSVResultParser(deconMSnResultFile);
		deconMSnResParser.parse(0.01f);
		Set<String> v2Scans = deconMSnResParser.getScanSet();
		
		System.out.println("ScanNum\tPrecursorMz(v1)\tCharge(v1)\tSpecEValue(v1)\tPrecursorMz(v2)\tCharge(v2)");
		
		for(String scan1 : msConvertScans)
		{
			if(!v2Scans.contains(scan1))
			{
				int scan = Integer.parseInt(scan1);
				String id1 = scanToId1.get(scan);
				Spectrum spec1 = v1SpecAccessor.getSpectrumById(id1);
				
				String id2 = scanToId2.get(scan);
				Spectrum spec2 = null;
				if(id2 != null)
					spec2 = deconMSnSpecAccessor.getSpectrumById(id2);
				
				if(spec2 != null)
				{
					System.out.println(spec1.getScanNum()
							+"\t"+spec1.getPrecursorPeak().getMz() + "\t" + spec1.getCharge() + "\t" + msConvertResParser.getSpecEValue(id1)
							+"\t"+spec2.getPrecursorPeak().getMz() + "\t" + spec2.getCharge()
							);
				}
				else
				{
					System.out.println(spec1.getScanNum()
							+"\t"+spec1.getPrecursorPeak().getMz() + "\t" + spec1.getCharge() + "\t" + msConvertResParser.getSpecEValue(id1)
							+"\t"+ "0" + "\t" + "0"
							);
				}
			}
		}
	}	
	
	@Test
	public void compareDeconMSnV1V3()
	{
		File dir = new File("C:\\cygwin\\home\\kims336\\Data\\IPA");

		File v1SpecFile = new File(dir.getPath() + File.separator + "QC_Shew_08_04-pt5-2_11jan09_sphinx_08-11-18_0515.mgf");
		SpectraAccessor v1SpecAccessor = new SpectraAccessor(v1SpecFile);
		
		File v1ResultFile = new File(dir.getPath() + File.separator + "Ver0515_10ppm_TI0_OnePerScan.tsv");
		TSVResultParser v1ResParser = new TSVResultParser(v1ResultFile);
		v1ResParser.parse(0.01f);
		Set<String> v1Scans = v1ResParser.getScanSet();
		Map<String,Integer> idToScan1 = new HashMap<String,Integer>();
		Map<Integer,String> scanToId1 = new HashMap<Integer,String>();
		Iterator<Spectrum> v1SpecItr = v1SpecAccessor.getSpecItr();
		while(v1SpecItr.hasNext())
		{
			Spectrum spec = v1SpecItr.next();
			idToScan1.put(spec.getID(), spec.getScanNum());
			scanToId1.put(spec.getScanNum(), spec.getID());
		}
		
		File v2SpecFile = new File(dir.getPath() + File.separator + "DeconMSn0516\\QC_Shew_08_04-pt5-2_11Jan09_Sphinx_08-11-18.mgf");
		SpectraAccessor v2SpecAccessor = new SpectraAccessor(v2SpecFile);

		Map<String,Integer> idToScan2 = new HashMap<String,Integer>();
		Map<Integer,String> scanToId2 = new HashMap<Integer,String>();
		Iterator<Spectrum> v2SpecItr = v2SpecAccessor.getSpecItr();
		while(v2SpecItr.hasNext())
		{
			Spectrum spec = v2SpecItr.next();
			idToScan2.put(spec.getID(), spec.getScanNum());
			scanToId2.put(spec.getScanNum(), spec.getID());
		}
		
		File v2ResultFile = new File(dir.getPath() + File.separator + "DeconMSn0516\\QC_Shew_08_04-pt5-2_11Jan09_Sphinx_08-11-18_OnePerScan.tsv");
		TSVResultParser v2ResParser = new TSVResultParser(v2ResultFile);
		v2ResParser.parse(0.01f);
		Set<String> v2Scans = v2ResParser.getScanSet();
		
		System.out.println("ScanNum\tPrecursorMz(v1)\tCharge(v1)\tSpecEValue(v1)\tPrecursorMz(v2)\tCharge(v2)");
		
		for(String scan1 : v1Scans)
		{
			if(!v2Scans.contains(scan1))
			{
				int scan = Integer.parseInt(scan1);
				String id1 = scanToId1.get(scan);
				Spectrum spec1 = v1SpecAccessor.getSpectrumById(id1);
				
				String id2 = scanToId2.get(scan);
				Spectrum spec2 = null;
				if(id2 != null)
					spec2 = v2SpecAccessor.getSpectrumById(id2);
				
				if(spec2 != null)
				{
					System.out.println(spec1.getScanNum()
							+"\t"+spec1.getPrecursorPeak().getMz() + "\t" + spec1.getCharge() + "\t" + v1ResParser.getSpecEValue(id1)
							+"\t"+spec2.getPrecursorPeak().getMz() + "\t" + spec2.getCharge()
							);
				}
				else
				{
					System.out.println(spec1.getScanNum()
							+"\t"+spec1.getPrecursorPeak().getMz() + "\t" + spec1.getCharge() + "\t" + v1ResParser.getSpecEValue(id1)
							+"\t"+ "0" + "\t" + "0"
							);
				}
			}
		}
	}
}
