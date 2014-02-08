package msgfplus;

import java.io.File;

import org.junit.Test;

import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.mzml.MzMLAdapter;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;
import edu.ucsd.msjava.ui.ScoringParamGen;

public class TestScoring {
	@Test
	public void testReadingParamFile()
	{
		String paramFile = "D:\\Research\\Data\\Ansong_TMT_Velos\\params\\HCD_HighRes_Tryp_TMT.param";
		NewRankScorer scorer = new NewRankScorer(paramFile);
	}
	
	@Test
	public void testScoringParamGen()
	{
		File resultPath = new File("D:\\Research\\Data\\Ansong_TMT_Velos\\TMT_QE");
		File specPath = new File("D:\\Research\\Data\\Ansong_TMT_Velos\\dta");

		String[] argv = {"-i", resultPath.getPath(), "-d", specPath.getPath(), "-m", "2", "-inst", "1", "-e", "0" 
				,"-protocol", "4"
				};
		
		ParamManager paramManager = new ParamManager("ScoringParamGen", "Test", "Test",
				"java -Xmx2000M -cp MSGFPlus.jar edu.ucsd.msjava.ui.ScoringParamGen");
			
		MzMLAdapter.turnOffLogs();
		paramManager.addScoringParamGenParams();
		paramManager.parseParams(argv);
		ScoringParamGen.runScoringParamGen(paramManager);
		System.out.println("Done");		
	}		
}
