package msdbsearch;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.TreeMap;

import sequences.Constants;

import msgf.NominalMass;
import msgf.Tolerance;
import msscorer.DBScanScorer;
import msscorer.FastScorer;
import msscorer.NewRankScorer;
import msscorer.NewScoredSpectrum;
import msscorer.NewScorerFactory;
import msutil.ActivationMethod;
import msutil.Composition;
import msutil.Enzyme;
import msutil.SpecKey;
import msutil.Spectrum;
import msutil.SpectrumAccessorByScanNum;

public class SpectraScanner {
	private TreeMap<Float,SpecKey> pepMassSpecKeyMap;
	private HashMap<SpecKey,FastScorer> specKeyScorerMap;

	public SpectraScanner(
			SpectrumAccessorByScanNum specMap,
			List<SpecKey> specKeyList,
    		Tolerance leftParentMassTolerance, 
    		Tolerance rightParentMassTolerance, 
			int numAllowedC13,
			NewRankScorer scorer,
			ActivationMethod activationMethod,
			Enzyme enzyme
			)
	{
		for(SpecKey specKey : specKeyList)
		{
			int scanNum = specKey.getScanNum();
			Spectrum spec = specMap.getSpectrumByScanNum(scanNum);
			if(activationMethod != null && spec.getActivationMethod() != null && (spec.getActivationMethod() != activationMethod))
				continue;
			if(spec.size() < Constants.MIN_NUM_PEAKS_PER_SPECTRUM)
			{
				System.out.println("Spectrum " + spec.getScanNum() + " has too few peaks (#Peaks: " + spec.size()+"): ignored.");
				continue;
			}
			if(scorer == null)
			{
				scorer = NewScorerFactory.get(spec.getActivationMethod(), enzyme);
			}
			
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
}
