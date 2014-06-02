package msgfplus;

import java.io.File;

import org.junit.Test;

import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.ScoringParameterGeneratorWithErrors;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.Protocol;
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
		File resultPath = new File("C:\\cygwin\\home\\kims336\\Data\\Scoring");
		File specPath = new File("C:\\cygwin\\home\\kims336\\Data\\Scoring");

		String[] argv = {"-i", resultPath.getPath(), "-d", specPath.getPath(), "-m", "2", "-inst", "3", "-e", "0" 
				,"-protocol", "5"
				};
		
		ParamManager paramManager = new ParamManager("ScoringParamGen", "Test", "Test",
				"java -Xmx2000M -cp MSGFPlus.jar edu.ucsd.msjava.ui.ScoringParamGen");
			
		MzMLAdapter.turnOffLogs();
		paramManager.addScoringParamGenParams();
		paramManager.parseParams(argv);
		ScoringParamGen.runScoringParamGen(paramManager);
		System.out.println("Done");		
	}		
	
	@Test
	public void testScoringParamGenFromMgf()
	{
		ActivationMethod actMethod = ActivationMethod.HCD;
		InstrumentType instType = InstrumentType.QEXACTIVE;
		Enzyme enzyme = Enzyme.TRYPSIN;
		Protocol protocol = Protocol.STANDARD;
		File specFile = new File("D:\\Research\\Data\\TrainingMSGFPlus\\AnnotatedSpectra\\HCD_QExactive_Tryp.mgf");
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		
		SpecDataType dataType = new SpecDataType(actMethod, instType, enzyme, protocol);
		System.out.println("Processing " + dataType.toString());
		ScoringParameterGeneratorWithErrors.generateParameters(
				specFile,
				dataType,
				aaSet, 
				new File("C:\\cygwin\\home\\kims336\\Data\\Scoring"),
				false, 
				false,
				false);
		
	}
}
