package msdbsearch;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import msgf.NominalMass;
import msgf.ScoredSpectrum;
import msgf.ScoredSpectrumSum;
import msgf.Tolerance;
import msscorer.DBScanScorer;
import msscorer.DBScanScorerSum;
import msscorer.FastScorer;
import msscorer.NewRankScorer;
import msscorer.NewScoredSpectrum;
import msscorer.NewScorerFactory;
import msscorer.SimpleDBSearchScorer;
import msscorer.NewScorerFactory.SpecDataType;
import msutil.ActivationMethod;
import msutil.Composition;
import msutil.Enzyme;
import msutil.InstrumentType;
import msutil.Modification;
import msutil.Pair;
import msutil.Protocol;
import msutil.SpecKey;
import msutil.Spectrum;
import msutil.SpectrumAccessorBySpecIndex;

public class ScoredSpectraMap {
	private final SpectrumAccessorBySpecIndex specMap;
	private final List<SpecKey> specKeyList;
	private final Tolerance leftParentMassTolerance;
	private final Tolerance rightParentMassTolerance;
	private final int minNum13C;
	private final int maxNum13C;
	private final SpecDataType specDataType;
	
	private SortedMap<Double,SpecKey> pepMassSpecKeyMap;
	private Map<SpecKey,SimpleDBSearchScorer<NominalMass>> specKeyScorerMap;
	private Map<Pair<Integer,Integer>, SpecKey> specIndexChargeToSpecKeyMap;
	
	private boolean turnOffEdgeScoring = false;

	public ScoredSpectraMap(
			SpectrumAccessorBySpecIndex specMap,
			List<SpecKey> specKeyList,			
    		Tolerance leftParentMassTolerance, 
    		Tolerance rightParentMassTolerance, 
			int minNum13C,
			int maxNum13C,
			SpecDataType specDataType
			)
	{
		this.specMap = specMap;
		this.specKeyList = specKeyList;
		this.leftParentMassTolerance = leftParentMassTolerance;
		this.rightParentMassTolerance = rightParentMassTolerance;
		this.minNum13C = minNum13C;
		this.maxNum13C = maxNum13C;
		this.specDataType = specDataType;
		
		pepMassSpecKeyMap = Collections.synchronizedSortedMap((new TreeMap<Double,SpecKey>()));
		specKeyScorerMap = Collections.synchronizedMap(new HashMap<SpecKey,SimpleDBSearchScorer<NominalMass>>());
		specIndexChargeToSpecKeyMap = Collections.synchronizedMap(new HashMap<Pair<Integer,Integer>,SpecKey>());
	}

	public ScoredSpectraMap(
			SpectrumAccessorBySpecIndex specMap,
			List<SpecKey> specKeyList,			
    		Tolerance leftParentMassTolerance, 
    		Tolerance rightParentMassTolerance, 
			int maxNum13C,
			SpecDataType specDataType
			)
	{
		this.specMap = specMap;
		this.specKeyList = specKeyList;
		this.leftParentMassTolerance = leftParentMassTolerance;
		this.rightParentMassTolerance = rightParentMassTolerance;
		this.minNum13C = 0;
		this.maxNum13C = maxNum13C;
		this.specDataType = specDataType;
		
		pepMassSpecKeyMap = Collections.synchronizedSortedMap((new TreeMap<Double,SpecKey>()));
		specKeyScorerMap = Collections.synchronizedMap(new HashMap<SpecKey,SimpleDBSearchScorer<NominalMass>>());
		specIndexChargeToSpecKeyMap = Collections.synchronizedMap(new HashMap<Pair<Integer,Integer>,SpecKey>());
	}
	
	public ScoredSpectraMap turnOffEdgeScoring()
	{
		this.turnOffEdgeScoring = true;
		return this;
	}
	
	public SortedMap<Double,SpecKey> getPepMassSpecKeyMap()		{ return pepMassSpecKeyMap; }
	public Map<SpecKey,SimpleDBSearchScorer<NominalMass>> getSpecKeyScorerMap()	{ return specKeyScorerMap; }
	public Tolerance getLeftParentMassTolerance()				{ return leftParentMassTolerance; }
	public Tolerance getRightParentMassTolerance()				{ return rightParentMassTolerance; }
//	public int getNumAllowedC13()								{ return numAllowedC13; }
	public int getMaxNum13C()										{ return maxNum13C; }
	public int getMinNum13C()										{ return minNum13C; }
	
	public List<SpecKey> getSpecKeyList()	{ return specKeyList; }
	public SpecKey getSpecKey(int specIndex, int charge)
	{
		return specIndexChargeToSpecKeyMap.get(new Pair<Integer,Integer>(specIndex, charge));
	}
	
	public ScoredSpectraMap makePepMassSpecKeyMap()
	{
		for(SpecKey specKey : specKeyList)
		{
			Spectrum spec = specMap.getSpectrumBySpecIndex(specKey.getSpecIndex());
			spec.setCharge(specKey.getCharge());
			float peptideMass = spec.getParentMass() - (float)Composition.H2O;
			double peptideMassKey = (double)peptideMass;
			while(pepMassSpecKeyMap.get(peptideMassKey) != null)	// for speeding up
				peptideMassKey = Math.nextUp(peptideMassKey);
			pepMassSpecKeyMap.put(peptideMassKey, specKey);
			
			for(int delta = this.minNum13C; delta<=maxNum13C; delta++)
			{
				float mass1 = peptideMass-delta*(float)Composition.ISOTOPE;
				double mass1Key = (double)mass1;
				while(pepMassSpecKeyMap.get(mass1Key) != null)
					mass1Key = Math.nextUp(mass1Key);
				pepMassSpecKeyMap.put(mass1Key, specKey);
			}
			specIndexChargeToSpecKeyMap.put(new Pair<Integer,Integer>(specKey.getSpecIndex(), specKey.getCharge()), specKey);
		}
		return this;
	}
	
	public void preProcessSpectra()
	{
		preProcessSpectra(0, specKeyList.size());
	}
	
	public void preProcessSpectra(int fromIndex, int toIndex)
	{
		if(specDataType.getActivationMethod() != ActivationMethod.FUSION)
			preProcessIndividualSpectra(fromIndex, toIndex);
		else
			preProcessFusedSpectra(fromIndex, toIndex);
	}
	
	private void preProcessIndividualSpectra(int fromIndex, int toIndex)
	{
		NewRankScorer scorer = null;
		ActivationMethod activationMethod = specDataType.getActivationMethod();
		InstrumentType instType = specDataType.getInstrumentType();
		Enzyme enzyme = specDataType.getEnzyme();
		Protocol protocol = specDataType.getProtocol();
		
		if(activationMethod != ActivationMethod.ASWRITTEN && activationMethod != ActivationMethod.FUSION)
		{
			scorer = NewScorerFactory.get(activationMethod, instType, enzyme, protocol);
			if(this.turnOffEdgeScoring)
				scorer.doNotUseError();
		}
		
		for(SpecKey specKey : specKeyList.subList(fromIndex, toIndex))
		{
			int specIndex = specKey.getSpecIndex();
			Spectrum spec = specMap.getSpectrumBySpecIndex(specIndex);
			if(activationMethod == ActivationMethod.ASWRITTEN || activationMethod == ActivationMethod.FUSION)
			{
				scorer = NewScorerFactory.get(spec.getActivationMethod(), instType, enzyme, protocol);
				if(this.turnOffEdgeScoring)
					scorer.doNotUseError();
			}
			int charge = specKey.getCharge();
			spec.setCharge(charge);
			NewScoredSpectrum<NominalMass> scoredSpec = scorer.getScoredSpectrum(spec);
			
			float peptideMass = spec.getParentMass() - (float)Composition.H2O;
			float tolDaLeft = leftParentMassTolerance.getToleranceAsDa(peptideMass);
			int maxNominalPeptideMass = NominalMass.toNominalMass(peptideMass) + Math.round(tolDaLeft-0.4999f) + 1;
			
			if(scorer.supportEdgeScores())
				specKeyScorerMap.put(specKey, new DBScanScorer(scoredSpec, maxNominalPeptideMass));
			else
				specKeyScorerMap.put(specKey, new FastScorer(scoredSpec, maxNominalPeptideMass));
		}				
	}
	
	private void preProcessFusedSpectra(int fromIndex, int toIndex)
	{
		InstrumentType instType = specDataType.getInstrumentType();
		Enzyme enzyme = specDataType.getEnzyme();
		Protocol protocol = specDataType.getProtocol();
		
		for(SpecKey specKey : specKeyList.subList(fromIndex, toIndex))
		{
			ArrayList<Integer> specIndexList = specKey.getSpecIndexList();
			if(specIndexList == null)
			{
				specIndexList = new ArrayList<Integer>();
				specIndexList.add(specKey.getSpecIndex());
			}
			ArrayList<ScoredSpectrum<NominalMass>> scoredSpecList = new ArrayList<ScoredSpectrum<NominalMass>>();
			boolean supportEdgeScore = true;
			for(int specIndex : specIndexList)
			{
				Spectrum spec = specMap.getSpectrumBySpecIndex(specIndex);
				
				NewRankScorer scorer = NewScorerFactory.get(spec.getActivationMethod(), instType, enzyme, protocol);
				if(!scorer.supportEdgeScores())
					supportEdgeScore = false;
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
			int maxNominalPeptideMass = NominalMass.toNominalMass(peptideMass) + Math.round(tolDaLeft-0.4999f) + 1;
			if(supportEdgeScore)
//				specKeyScorerMap.put(specKey, new DBScanScorerSum(scoredSpecList, maxNominalPeptideMass));
				specKeyScorerMap.put(specKey, new FastScorer(scoredSpec, maxNominalPeptideMass));
			else
				specKeyScorerMap.put(specKey, new FastScorer(scoredSpec, maxNominalPeptideMass));
		}				
	}
}
