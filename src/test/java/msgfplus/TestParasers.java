package msgfplus;

import org.junit.Test;

import edu.ucsd.msjava.msutil.SpectraMap;
import edu.ucsd.msjava.parser.MS2SpectrumParser;

public class TestParasers {
	@Test
	public void testMS2Parser()
	{
		String fileName = System.getProperty("user.home")+"/Research/Data/Viktor/QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.ms2";
		
		SpectraMap map = new SpectraMap(fileName, new MS2SpectrumParser());
		map.getSpectrumBySpecIndex(0);
		
	}
}
