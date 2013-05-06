package msgfplus;

import java.io.File;
import java.io.FileNotFoundException;

import org.junit.Test;

import edu.ucsd.msjava.misc.MS2ToMgf;
import edu.ucsd.msjava.misc.ParamToTxt;
import edu.ucsd.msjava.msutil.SpectraIterator;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.parser.MgfSpectrumParser;
import edu.ucsd.msjava.parser.PNNLSpectraIterator;
import edu.ucsd.msjava.parser.PNNLSpectrumParser;

public class TestParsers {
	@Test
	public void testParamToTxt()
	{
//		File paramFile = new File("/Users/kims336/Research/Data/TrainingITRAQ/Global/params/HCD_HighRes_Tryp_iTRAQ.param");
//		ParamToTxt.convert(paramFile);
	}
	
	@Test
	public void ms2ToMgfTest()
	{
//		File ms2File = new File(System.getProperty("user.home")+"/Research/Data/Viktor/103111-Yeast-2hr-01.ANNOTATED.ms2");
//		File mgfFile = new File(System.getProperty("user.home")+"/Research/Data/Viktor/103111-Yeast-2hr-01.ANNOTATED.mgf");
//		try {
//			MS2ToMgf.convert(ms2File, mgfFile);
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
	}
		
	@Test
	public void newDeconMSnTest()
	{
		int numSpecs = 0;
		File mgfFile = new File("C:\\cygwin\\home\\kims336\\Data\\IPA\\QC_Shew_08_04-pt5-2_11Jan09_Sphinx_08-11-18_version2.mgf");
		try {
			SpectraIterator itr = new SpectraIterator(mgfFile.getPath(), new MgfSpectrumParser());
			while(itr.hasNext())
			{
				Spectrum spec = itr.next();
				float precursorMass = spec.getPrecursorPeak().getMass();
				if(precursorMass <= 0)
				{
					System.out.println("ScanNum " + spec.getScanNum() + " " + precursorMass);
				}
				if(spec.getPrecursorTolerance() == null)
				{
					System.out.println("ScanNum " + spec.getScanNum() + " tolerance null.");
				}
				++numSpecs;
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		System.out.println("NumSpectra: " + numSpecs);
	}
	
	@Test
	public void testReadingDtaTxt()
	{
		int numSpecs = 0;
		File mgfFile = new File("C:\\cygwin\\home\\kims336\\Data\\QCShewQE\\QC_Shew_13_02_2500ng_B_18Mar13_Jaguar_13-03-11_dta.txt");
		try {
			SpectraIterator itr = new PNNLSpectraIterator(mgfFile.getPath());
			while(itr.hasNext())
			{
				Spectrum spec = itr.next();
				boolean isHighPrecision = spec.isHighPrecision();
				if(!isHighPrecision)
				{
					System.out.println("ScanNum " + spec.getScanNum());
				}
				++numSpecs;
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		System.out.println("NumSpectra: " + numSpecs);
	}
	
}
