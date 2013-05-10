package msgfplus;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.parser.TSVResultParser;
import edu.ucsd.msjava.ui.MSGFPlus;

public class TestIPA {
	@Test
	public void iPRG2013()
	{
		File dir = new File("D:\\Research\\Data\\iPRG2013");

		File specFile = new File(dir.getPath()+File.separator+"DeconMSn\\f01.mgf");
		File dbFile = new File(dir.getPath()+File.separator+"Homo_sapiens_non-redundant.GRCh37.68.pep.all_FPKM_NOVEL-cRAP_targetdecoy.fasta");
		File modFile = new File(dir.getPath()+File.separator+"Mods.txt");

		String[] argv = {"-s", specFile.getPath(), "-d", dbFile.getPath(), 
				"-mod", modFile.getPath(), "-t", "10ppm", "-tda", "0", "-ti", "0,0", "-inst" , "1", "-protocol", "4",
		}; 

		ParamManager paramManager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
		paramManager.addMSGFPlusParams();
		
		String msg = paramManager.parseParams(argv);
		assertTrue(msg == null);
		
		assertTrue(MSGFPlus.runMSGFPlus(paramManager) == null);
		
	}
	
	@Test
	public void compareDeconMSnV1V3()
	{
		File dir = new File("C:\\cygwin\\home\\kims336\\Data\\IPA");

		File v1SpecFile = new File(dir.getPath() + File.separator + "QC_Shew_08_04-pt5-2_11Jan09_Sphinx_08-11-18_version1_dta.txt");
		SpectraAccessor v1SpecAccessor = new SpectraAccessor(v1SpecFile);
		
		File v1ResultFile = new File(dir.getPath() + File.separator + "Ver1_10ppm_TI0.tsv");
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
		
		File v2SpecFile = new File(dir.getPath() + File.separator + "QC_Shew_08_04-pt5-2_11Jan09_Sphinx_08-11-18_version3.mgf");
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
		
		File v2ResultFile = new File(dir.getPath() + File.separator + "Ver3_10ppm_TI0.tsv");
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
