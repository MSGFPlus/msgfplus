package msgfplus;

import java.io.File;
import java.util.HashSet;
import java.util.Iterator;

import org.junit.Test;

import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.IonProbability;
import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.ScoringParameterGeneratorWithErrors;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.IonType;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Peptide;
import edu.ucsd.msjava.msutil.Protocol;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.TopNFilter;
import edu.ucsd.msjava.msutil.WindowFilter;
import edu.ucsd.msjava.msutil.IonType.PrefixIon;
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
	
	@Test
	public void testGeneratingIonProbabilities()
	{
		File mgfFile = new File("D:\\Research\\Data\\TrainingMSGFPlus\\AnnotatedSpectra\\CID_LowRes_Tryp.mgf");
		SpectraAccessor accessor = new SpectraAccessor(mgfFile);
		IonType[] prefixIons = {IonType.getIonType("b"), IonType.getIonType("b-H2O"), IonType.getIonType("b-NH3"),
				IonType.getIonType("a"),
				IonType.getIonType("b2"), IonType.getIonType("b2-H2O"), IonType.getIonType("b2-NH3"),
				IonType.getIonType("b3"), IonType.getIonType("b3-H2O"), IonType.getIonType("b3-NH3")};

		IonType[] suffixIons = {IonType.getIonType("y"), IonType.getIonType("y-H2O"), IonType.getIonType("y-NH3"),
				IonType.getIonType("y2"), IonType.getIonType("y2-H2O"), IonType.getIonType("y2-NH3"),
				IonType.getIonType("y3"), IonType.getIonType("y3-H2O"), IonType.getIonType("y3-NH3")};
		
		TopNFilter filter = new TopNFilter(10);

		HashSet<String> pepSet = new HashSet<String>();

		int numPrefix = 0;
		int numSuffix = 0;
		int numNoise = 0;
		int numSpectra = 0;
		
		int sumNumPeaks = 0;
		
		Iterator<Spectrum> itr = accessor.getSpecItr();
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			if(filter != null)
				spec = filter.apply(spec);
			
			Peptide pep = spec.getAnnotation();
			if(pep == null)
				continue;
			
			String pepStr = spec.getAnnotationStr();
			if(pepSet.contains(pepStr))
				continue;
			else
				pepSet.add(pepStr);
			
			++numSpectra;
			sumNumPeaks += spec.size();
			for(Peak p : spec)
			{
				double peakMz = p.getMz();
				
				boolean detected = false;
				
				// suffix
				double srm = 0;
				for(int i=0; i<pep.size(); i++)
				{
					srm += pep.get(pep.size()-1-i).getMass();
					for(IonType ion : suffixIons)
					{
						float mz = ion.getMz((float)srm);
						if(mz > peakMz - 0.5 && mz < peakMz + 0.5)
						{
							++numSuffix;
							detected = true;
							break;
						}
					}
					if(detected) break;
				}
				
				if(detected) continue;
				
				// prefix
				double prm = 0;
				for(int i=0; i<pep.size(); i++)
				{
					prm += pep.get(i).getMass();
					for(IonType ion : prefixIons)
					{
						float mz = ion.getMz((float)prm);
						if(mz > peakMz - 0.5 && mz < peakMz + 0.5)
						{
							++numPrefix;
							detected = true;
							break;
						}
					}
					if(detected) break;
				}
				if(detected) continue;
				
				++numNoise;
			}
		}
			
		int numPeaks = numPrefix + numSuffix + numNoise;
		System.out.println("NumSpectra: " + numSpectra);
		System.out.println("NumSpecPeaks: " + sumNumPeaks);
		System.out.println("NumPeaks: " + numPeaks);
		System.out.println("% Prefix: " + numPrefix/(float)numPeaks);
		System.out.println("% Suffix: " + numSuffix/(float)numPeaks);
		System.out.println("% Noise: " + numNoise/(float)numPeaks);
	}
}
