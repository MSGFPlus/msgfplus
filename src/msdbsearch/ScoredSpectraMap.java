package msdbsearch;

import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.TreeMap;

import sequences.Constants;

import msgf.NominalMass;
import msgf.ScoredSpectrum;
import msgf.ScoredSpectrumSum;
import msgf.Tolerance;
import msscorer.DBScanScorer;
import msscorer.FastScorer;
import msscorer.NewRankScorer;
import msscorer.NewScoredSpectrum;
import msscorer.NewScorerFactory;
import msscorer.SimpleDBSearchScorer;
import msutil.ActivationMethod;
import msutil.Composition;
import msutil.Enzyme;
import msutil.InstrumentType;
import msutil.SpecKey;
import msutil.Spectrum;
import msutil.SpectrumAccessorByScanNum;

public class ScoredSpectraMap {
	private Tolerance leftParentMassTolerance;
	private Tolerance rightParentMassTolerance;
	private int numAllowedC13;

	private TreeMap<Float,SpecKey> pepMassSpecKeyMap;
	private HashMap<SpecKey,SimpleDBSearchScorer<NominalMass>> specKeyScorerMap;

	public ScoredSpectraMap(
			SpectrumAccessorByScanNum specMap,
			List<SpecKey> specKeyList,
    		Tolerance leftParentMassTolerance, 
    		Tolerance rightParentMassTolerance, 
			int numAllowedC13,
			ActivationMethod activationMethod,
			InstrumentType instType,
			Enzyme enzyme
			)
	{
		this.leftParentMassTolerance = leftParentMassTolerance;
		this.rightParentMassTolerance = rightParentMassTolerance;
		this.numAllowedC13 = numAllowedC13;
		
		pepMassSpecKeyMap = new TreeMap<Float,SpecKey>();
		specKeyScorerMap = new HashMap<SpecKey,SimpleDBSearchScorer<NominalMass>>();
		
		if(activationMethod != ActivationMethod.FUSION)
			preProcessSpectra(specMap, specKeyList, activationMethod, instType, enzyme);
		else
			preProcessFusedSpectra(specMap, specKeyList, instType, enzyme);
	}
	
	public void preProcessSpectra(
			SpectrumAccessorByScanNum specMap,
			List<SpecKey> specKeyList,
			ActivationMethod activationMethod,
			InstrumentType instType,
			Enzyme enzyme
			)
	{
		NewRankScorer scorer = null;
		if(activationMethod != null && activationMethod != ActivationMethod.FUSION)
			scorer = NewScorerFactory.get(activationMethod, instType, enzyme);
		
		for(SpecKey specKey : specKeyList)
		{
			int scanNum = specKey.getScanNum();
			Spectrum spec = specMap.getSpectrumByScanNum(scanNum);
			if(spec.size() < Constants.MIN_NUM_PEAKS_PER_SPECTRUM)
			{
//				System.out.println("Spectrum " + spec.getScanNum() + " has too few peaks (#Peaks: " + spec.size()+"): ignored.");
				continue;
			}
			if(activationMethod == null || activationMethod == ActivationMethod.FUSION)
				scorer = NewScorerFactory.get(spec.getActivationMethod(), instType, enzyme);
			
			int charge = specKey.getCharge();
			spec.setCharge(charge);
			NewScoredSpectrum<NominalMass> scoredSpec = scorer.getScoredSpectrum(spec);
			
			float peptideMass = spec.getParentMass() - (float)Composition.H2O;
			float tolDaLeft = leftParentMassTolerance.getToleranceAsDa(peptideMass);
			int maxNominalPeptideMass = NominalMass.toNominalMass(peptideMass) + Math.round(tolDaLeft-0.4999f);
			
			if(scorer.supportEdgeScores())
				specKeyScorerMap.put(specKey, new DBScanScorer(scoredSpec, maxNominalPeptideMass));
			else
				specKeyScorerMap.put(specKey, new FastScorer(scoredSpec, maxNominalPeptideMass));
			
			while(pepMassSpecKeyMap.get(peptideMass) != null)	// for speeding up
				peptideMass = Math.nextUp(peptideMass);
			pepMassSpecKeyMap.put(peptideMass, specKey);
			
			float tolDaRight = rightParentMassTolerance.getToleranceAsDa(peptideMass);
			if(numAllowedC13 > 0 && tolDaRight < 0.5f)
			{
				if(numAllowedC13 >= 1)
				{
					float mass1 = peptideMass-(float)Composition.ISOTOPE;
					while(pepMassSpecKeyMap.get(mass1) != null)
						mass1 = Math.nextUp(mass1);
					pepMassSpecKeyMap.put(mass1, specKey);
				}
				
				if(numAllowedC13 >= 2)
				{
					float mass2 = peptideMass-2*(float)Composition.ISOTOPE;
					while(pepMassSpecKeyMap.get(mass2) != null)
						mass2 = Math.nextUp(mass2);
					pepMassSpecKeyMap.put(mass2, specKey);
				}
			}				
		}				
	}
	
	public void preProcessFusedSpectra(
			SpectrumAccessorByScanNum specMap,
			List<SpecKey> specKeyList,
			InstrumentType instType,
			Enzyme enzyme
			)
	{
		for(SpecKey specKey : specKeyList)
		{
			ArrayList<Integer> scanNumList = specKey.getScanNumList();
			if(scanNumList == null)
			{
				scanNumList = new ArrayList<Integer>();
				scanNumList.add(specKey.getScanNum());
			}
			ArrayList<ScoredSpectrum<NominalMass>> scoredSpecList = new ArrayList<ScoredSpectrum<NominalMass>>();
			for(int scanNum : scanNumList)
			{
				Spectrum spec = specMap.getSpectrumByScanNum(scanNum);
				if(spec.size() < Constants.MIN_NUM_PEAKS_PER_SPECTRUM)
				{
//					System.out.println("Spectrum " + spec.getScanNum() + " has too few peaks (#Peaks: " + spec.size()+"): ignored.");
					continue;
				}
				
				NewRankScorer scorer = NewScorerFactory.get(spec.getActivationMethod(), instType, enzyme);
				int charge = specKey.getCharge();
				spec.setCharge(charge);
				NewScoredSpectrum<NominalMass> sSpec = scorer.getScoredSpectrum(spec);
				scoredSpecList.add(sSpec);
			}
			
			if(scoredSpecList.size() == 0)
				continue;
			ScoredSpectrumSum<NominalMass> scoredSpec = new ScoredSpectrumSum<NominalMass>(scoredSpecList);
			float peptideMass = scoredSpec.getPrecursorPeak().getMass() - (float)Composition.H2O;
			float tolDaLeft = leftParentMassTolerance.getToleranceAsDa(peptideMass);
			int maxNominalPeptideMass = NominalMass.toNominalMass(peptideMass) + Math.round(tolDaLeft-0.4999f);
			specKeyScorerMap.put(specKey, new FastScorer(scoredSpec, maxNominalPeptideMass));
			
			while(pepMassSpecKeyMap.get(peptideMass) != null)	// for speeding up
				peptideMass = Math.nextUp(peptideMass);
			pepMassSpecKeyMap.put(peptideMass, specKey);
			
			float tolDaRight = rightParentMassTolerance.getToleranceAsDa(peptideMass);
			if(numAllowedC13 > 0 && tolDaRight < 0.5f)
			{
				if(numAllowedC13 >= 1)
				{
					float mass1 = peptideMass-(float)Composition.ISOTOPE;
					while(pepMassSpecKeyMap.get(mass1) != null)
						mass1 = Math.nextUp(mass1);
					pepMassSpecKeyMap.put(mass1, specKey);
				}
				
				if(numAllowedC13 >= 2)
				{
					float mass2 = peptideMass-2*(float)Composition.ISOTOPE;
					while(pepMassSpecKeyMap.get(mass2) != null)
						mass2 = Math.nextUp(mass2);
					pepMassSpecKeyMap.put(mass2, specKey);
				}
			}					
		}				
	}
	
	public TreeMap<Float,SpecKey> getPepMassSpecKeyMap()		{ return pepMassSpecKeyMap; }
	public HashMap<SpecKey,SimpleDBSearchScorer<NominalMass>> getSpecKeyScorerMap()	{ return specKeyScorerMap; }
	public Tolerance getLeftParentMassTolerance()				{ return leftParentMassTolerance; }
	public Tolerance getRightParentMassTolerance()				{ return rightParentMassTolerance; }
	public int getNumAllowedC13()								{ return numAllowedC13; }
	
}
