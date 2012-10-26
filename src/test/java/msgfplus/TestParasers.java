package msgfplus;

import java.io.File;

import org.junit.Test;

import edu.ucsd.msjava.misc.ParamToTxt;

public class TestParasers {
	@Test
	public void testParamToTxt()
	{
		File paramFile = new File("/Users/kims336/Research/Data/TrainingITRAQ/Global/params/HCD_HighRes_Tryp_iTRAQ.param");
		ParamToTxt.convert(paramFile);
	}
}
