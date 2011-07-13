package msgf.test;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import msdbsearch.DBScanner;
import msgf.GeneratingFunction;
import msgf.GenericDeNovoGraph;
import msgf.IntMassFactory;
import msgf.NominalMass;
import msgf.NominalMassFactory;
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
import msutil.SpectrumAccessorBySpecIndex;

import parser.InsPecTPSM;
import parser.InsPecTParser;
import parser.MgfSpectrumParser;
import parser.PSMList;

public class MSGFTest {
	public static void main(String[] argv) throws Exception
	{
//		msgfTest();
		dicTest();
	}
	
	public static void dicTest() throws Exception
	{
		File specFile = new File("/home/sangtaekim/Research/Data/Zubarev/SACTest/SACTest.mgf");
		SpectrumAccessorBySpecIndex specAccessor = new SpectraMap(specFile.getPath(), new MgfSpectrumParser());
		int scanNum = 3888;
		Spectrum spec = specAccessor.getSpectrumBySpecIndex(scanNum);
		
		Enzyme enzyme = Enzyme.TRYPSIN;
//		AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile("/home/sangtaekim/Developments/MS_Java_Dev/bin/Mods.txt");
		AminoAcidSet aaSet = AminoAcidSet.getStandardAminoAcidSetWithFixedCarbamidomethylatedCys();
		DBScanner.setAminoAcidProbabilities("/home/sangtaekim/Research/Data/CommonContaminants/IPI_human_3.79_withContam.fasta", aaSet);
		aaSet.registerEnzyme(enzyme);
//		aaSet.printAASet();
		NewRankScorer scorer = NewScorerFactory.get(ActivationMethod.HCD, enzyme);
		
		Annotation annotation = new Annotation(getAnnotationStr("R.TPEVTCVVVDVSHEDPEVQFK.W"), aaSet);
		NewScoredSpectrum<NominalMass> scoredSpec1 = scorer.getScoredSpectrum(spec);
//		ScoredSpectrum<NominalMass> scoredSpec = new msscorer.DBScanScorer(scoredSpec1, annotation.getPeptide().getNominalMass());
		NominalMassFactory factory = new NominalMassFactory(aaSet, enzyme, 50);
		Tolerance pmTolerance = Tolerance.ZERO_TOLERANCE;
		GenericDeNovoGraph<NominalMass> graph = new GenericDeNovoGraph<NominalMass>(factory, spec.getParentMass(), pmTolerance, enzyme, scoredSpec1);
		GeneratingFunction<NominalMass> gf = new GeneratingFunction<NominalMass>(graph).enzyme(enzyme).doNotCalcNumber();
		gf.computeGeneratingFunction();
		int score = gf.getScore(annotation);
		double specProb = gf.getSpectralProbability(score);
		System.out.println(scanNum+"\t"+annotation+" "+(gf.getMaxScore()-1)+" "+score+" "+specProb);
		for(int t=score; t<gf.getMaxScore(); t++)
			System.out.println(t+"\t"+gf.getScoreDist().getProbability(t));
		
		System.out.println("SpecProb from Dictionary");
		ArrayList<String> dictionary = gf.getReconstructionsEqualOrAboveScore(score);
		float specProb2 = 0;
		HashMap<Integer,Float> scoreMap = new HashMap<Integer,Float>();
		for(int t=score; t<gf.getMaxScore(); t++)
			scoreMap.put(t, 0f);
		for(String pep : dictionary)
		{
			int t = gf.getScore(new Annotation(getAnnotationStr(pep), aaSet));
			float prob = getProbability(pep, aaSet);
			scoreMap.put(t, scoreMap.get(t)+prob);
			specProb2 += prob;
			System.out.println(pep+"\t"+prob+"\t"+t+"\t"+specProb2);
		}
		System.out.println("SpecProb2\t"+specProb2+"\t"+dictionary.size());
		for(int t=score; t<gf.getMaxScore(); t++)
			System.out.println(t+"\t"+scoreMap.get(t));
	}
	
	public static float getProbability(String annotationStr, AminoAcidSet aaSet)
	{
		float prob = 1;
		char aaBefore = annotationStr.charAt(0);
		if(aaBefore == 'K' || aaBefore == 'R' || !Character.isLetter(aaBefore))
			prob *= (aaSet.getAminoAcid('K').getProbability()+aaSet.getAminoAcid('R').getProbability());
		else
			prob *= (1-(aaSet.getAminoAcid('K').getProbability()+aaSet.getAminoAcid('R').getProbability()));
		
		String pepStr = annotationStr.substring(annotationStr.indexOf('.')+1);
		for(int i=0; i<pepStr.length(); i++)
		{
			char c = pepStr.charAt(i);
			if(Character.isLetter(c))
			{
				prob *= aaSet.getAminoAcid(Character.toUpperCase(c)).getProbability();
			}
		}
		return prob;
	}
	
	public static String getAnnotationStr(String annotationStr)
	{
		String retStr = null;
		retStr = annotationStr.replaceAll("E-18\\.011", "e");
		retStr = retStr.replaceAll("Q-17\\.027", "q");
		retStr = retStr.replaceAll("M\\+15\\.995", "m");
		if(annotationStr.charAt(annotationStr.length()-2) != '.')
			return retStr+".A";
		else
			return retStr;
	}
	
	
	public static void msgfTest() throws Exception
	{
		File specFile = new File(System.getProperty("user.home")+"/Research/Data/ISBControl/Mix_7/ORBITRAP/mgf/ISB02_mAB_ChymoTryp_Tryp.mgf");
		SpectrumAccessorBySpecIndex specAccessor = new SpectraMap(specFile.getPath(), new MgfSpectrumParser());
		
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
//			if(psm.getScanNum() != 1649)
//				continue;
			if(psm.getScanNum() != 7400)
			continue;
			Spectrum spec = specAccessor.getSpectrumBySpecIndex(psm.getScanNum());
			
			NewScoredSpectrum<IntMass> scoredSpec = scorer.getScoredSpectrum(spec);
			GenericDeNovoGraph<IntMass> graph;
			graph = new GenericDeNovoGraph<IntMass>(factory, spec.getParentMass(), pmTolerance, enzyme, scoredSpec);
			GeneratingFunction<IntMass> gf = new GeneratingFunction<IntMass>(graph).enzyme(enzyme);
			gf.computeGeneratingFunction();
			double specProb = gf.getSpectralProbability(psm.getAnnotation());
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
						prob *= aaSet.getProbCleavageSites();
					else
						prob *= 1-aaSet.getProbCleavageSites();
				}
				sumProb += prob;
				System.out.println(annotationStr+"\t"+prob+"\t"+sumProb+"\t"+gf.getScore(annotation)+"\t"+(pep.getMass()-spec.getParentMass()+(float)Composition.H2O));
			}
			
			System.out.println("SpecProb\t"+specProb+"\t"+sumProb);
//			System.exit(0);
		}
	}
}
