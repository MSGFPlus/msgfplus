package edu.ucsd.msjava.msdbsearch;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msgf.ScoredSpectrum;
import edu.ucsd.msjava.msgf.ScoredSpectrumSum;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msscorer.DBScanScorer;
import edu.ucsd.msjava.msscorer.FastScorer;
import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.NewScoredSpectrum;
import edu.ucsd.msjava.msscorer.NewScorerFactory;
import edu.ucsd.msjava.msscorer.SimpleDBSearchScorer;
import edu.ucsd.msjava.msscorer.NewScorerFactory.SpecDataType;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.Pair;
import edu.ucsd.msjava.msutil.Protocol;
import edu.ucsd.msjava.msutil.SpecKey;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;

public class ScoredSpectraMap {
	private final SpectraAccessor specAcc;
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
			SpectraAccessor specAcc,
			List<SpecKey> specKeyList,			
    		Tolerance leftParentMassTolerance, 
    		Tolerance rightParentMassTolerance, 
			int minNum13C,
			int maxNum13C,
			SpecDataType specDataType
			)
	{
		this.specAcc = specAcc;
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
			SpectraAccessor specAcc,
			List<SpecKey> specKeyList,			
    		Tolerance leftParentMassTolerance, 
    		Tolerance rightParentMassTolerance, 
			int maxNum13C,
			SpecDataType specDataType
			)
	{
		this(specAcc, specKeyList, leftParentMassTolerance, rightParentMassTolerance, 0, maxNum13C, specDataType);
	}
	
	public ScoredSpectraMap turnOffEdgeScoring()
	{
		this.turnOffEdgeScoring = true;
		return this;
	}
	
	public SortedMap<Double,SpecKey> getPepMassSpecKeyMap()		{ return pepMassSpecKeyMap; }
	public Map<SpecKey,SimpleDBSearchScorer<NominalMass>> getSpecKeyScorerMap()	{ return specKeyScorerMap; }
	public SpectraAccessor	getSpectraAccessor()				{ return specAcc; }
	public SpecDataType getSpecDataType()						{ return specDataType; }
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
			int specIndex = specKey.getSpecIndex();
			Spectrum spec = specAcc.getSpectrumBySpecIndex(specIndex);
			spec.setCharge(specKey.getCharge());
			float peptideMass = spec.getParentMass() - (float)Composition.H2O;
			
			for(int delta = this.minNum13C; delta<=maxNum13C; delta++)
			{
				float mass1 = peptideMass-delta*(float)Composition.ISOTOPE;
				double mass1Key = (double)mass1;
				while(pepMassSpecKeyMap.get(mass1Key) != null)
					mass1Key = Math.nextUp(mass1Key);
				pepMassSpecKeyMap.put(mass1Key, specKey);
			}
			specIndexChargeToSpecKeyMap.put(new Pair<Integer,Integer>(specIndex, specKey.getCharge()), specKey);
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
			Spectrum spec = specAcc.getSpectrumBySpecIndex(specIndex);
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
				Spectrum spec = specAcc.getSpectrumBySpecIndex(specIndex);
				
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
