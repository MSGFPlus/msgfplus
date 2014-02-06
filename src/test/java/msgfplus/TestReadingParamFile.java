package msgfplus;

import org.junit.Test;

import edu.ucsd.msjava.msscorer.NewRankScorer;

public class TestReadingParamFile {
	@Test
	public void testReadingParamFile()
	{
		String paramFile = "D:\\Research\\Data\\Ansong_TMT_Velos\\params\\HCD_HighRes_Tryp_TMT.param";
		NewRankScorer scorer = new NewRankScorer(paramFile);
	}
}
