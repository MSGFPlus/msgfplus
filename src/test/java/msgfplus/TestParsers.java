package msgfplus;

import java.io.File;

import org.junit.Test;

import edu.ucsd.msjava.misc.MS2ToMgf;
import edu.ucsd.msjava.misc.ParamToTxt;

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
		File ms2File = new File(System.getProperty("user.home")+"/Research/Data/Viktor/103111-Yeast-2hr-01.ANNOTATED.ms2");
		File mgfFile = new File(System.getProperty("user.home")+"/Research/Data/Viktor/103111-Yeast-2hr-01.ANNOTATED.mgf");
		try {
			MS2ToMgf.convert(ms2File, mgfFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
