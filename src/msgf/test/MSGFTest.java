package msgf.test;

import java.io.File;
import java.util.ArrayList;

import msgf.GeneratingFunction;
import msgf.GenericDeNovoGraph;
import msgf.IntMassFactory;
import msgf.ScoredSpectrum;
import msgf.Tolerance;
import msgf.IntMassFactory.IntMass;
import msscorer.NewRankScorer;
import msscorer.NewScoredSpectrum;
import msscorer.NewScorerFactory;
import msutil.ActivationMethod;
import msutil.AminoAcidSet;
import msutil.Annotation;
import msutil.Composition;
import msutil.Constants;
import msutil.Enzyme;
import msutil.Peptide;
import msutil.SpectraMap;
import msutil.Spectrum;
import msutil.SpectrumAccessorByScanNum;

import parser.InsPecTPSM;
import parser.InsPecTParser;
import parser.MgfSpectrumParser;
import parser.PSMList;

public class MSGFTest {
	public static void main(String[] argv) throws Exception
	{
		msgfTest();
	}
	
	public static void msgfTest() throws Exception
	{
		File specFile = new File(System.getProperty("user.home")+"/Research/Data/ISBControl/Mix_7/ORBITRAP/mgf/ISB02_mAB_ChymoTryp_Tryp.mgf");
		SpectrumAccessorByScanNum specAccessor = new SpectraMap(specFile.getPath(), new MgfSpectrumParser());
		
		File resultFile = new File(System.getProperty("user.home")+"/Research/Data/TDATest/MSGFDBv2/MSGFDB_S10_Target_New.txt");
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		InsPecTParser parser = new InsPecTParser(aaSet);
		parser.parse(resultFile.getPath());
		
		int maxLength = 40;
		float rescalingFactor = Constants.INTEGER_MASS_SCALER;
		rescalingFactor = Constants.INTEGER_MASS_SCALER_HIGH_PRECISION;
		Tolerance pmTolerance = Tolerance.ZERO_TOLERANCE;
		pmTolerance = new Tolerance(30, true);
		Enzyme enzyme = Enzyme.TRYPSIN;
//		enzyme = null;
		IntMassFactory factory = new IntMassFactory(aaSet, enzyme, maxLength, rescalingFactor);
		NewRankScorer scorer = NewScorerFactory.get(ActivationMethod.CID, enzyme);
		
		String header = parser.getHeader();
		PSMList<InsPecTPSM> psmList = parser.getPSMList(); 
		for(InsPecTPSM psm : psmList)
		{
//			if(psm.getPeptide().size() > 10)
//				continue;
			if(psm.getScanNum() != 1649)
				continue;
			Spectrum spec = specAccessor.getSpectrumByScanNum(psm.getScanNum());
			
			NewScoredSpectrum<IntMass> scoredSpec = scorer.getScoredSpectrum(spec);
			GenericDeNovoGraph<IntMass> graph;
			graph = new GenericDeNovoGraph<IntMass>(factory, spec.getParentMass(), pmTolerance, enzyme, scoredSpec);
			GeneratingFunction<IntMass> gf = new GeneratingFunction<IntMass>(graph).enzyme(enzyme);
			gf.computeGeneratingFunction();
			float specProb = gf.getSpectralProbability(psm.getAnnotation());
			int score = gf.getScore(psm.getAnnotation());
//			if(score != gf.getMaxScore()-1)
//				continue;
			
			int msgfScore = gf.getScore(psm.getAnnotation());
			ArrayList<String> dictionary = gf.getReconstructionsEqualOrAboveScore(msgfScore);
			
			System.out.println(psm.getAnnotation()+"\t"+msgfScore+"\t"+(gf.getMaxScore()-1));
			float sumProb = 0;
			for(String annotationStr : dictionary)
			{
				if(enzyme == null)
					annotationStr = "."+annotationStr;
				Annotation annotation = new Annotation(annotationStr+".A", aaSet);
				Peptide pep = annotation.getPeptide();
//				Peptide pep = aaSet.getPeptide(annotationStr);
				float prob = pep.getProbability();
				if(enzyme != null)
				{
					if(annotation.toString().startsWith("R."))
						prob *= enzyme.getProbCleavageSites();
					else
						prob *= 1-enzyme.getProbCleavageSites();
				}
				sumProb += prob;
				System.out.println(annotationStr+"\t"+prob+"\t"+sumProb+"\t"+gf.getScore(annotation)+"\t"+(pep.getMass()-spec.getParentMass()+(float)Composition.H2O));
			}
			
			System.out.println("SpecProb\t"+specProb+"\t"+sumProb);
//			System.exit(0);
		}
	}
}
