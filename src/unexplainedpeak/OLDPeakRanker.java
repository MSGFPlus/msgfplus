package unexplainedpeak;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import msgf.Tolerance;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Composition;
import msutil.IonType;
import msutil.Peak;
import msutil.Peptide;
import msutil.Spectrum;

public class OLDPeakRanker {
	AminoAcidSet aaSet;
	static float maxAAMass = 0, minAAMass = Float.MAX_VALUE;
	static Tolerance tolerance = new Tolerance(0.5f, false);
	
	static ArrayList<IonType> precursorIons = null;
	
	static void registerPrecursorIons(ArrayList<IonType> sigIons){
		precursorIons = new  ArrayList<IonType>();
		
		for(IonType ion : sigIons){
			if(ion instanceof IonType.PrecursorIon)
				precursorIons.add(ion);
		}
	}
	
	
	public OLDPeakRanker(AminoAcidSet aaSet, Tolerance tol){
		this.aaSet = aaSet;
		tolerance = tol;
		if(maxAAMass == 0){
			for(AminoAcid aa : aaSet){
				minAAMass = Math.min(aa.getMass(), minAAMass);
				maxAAMass = Math.max(aa.getMass(), maxAAMass);
			}
		}
	}
	
	
	HashMap<Integer, Peak> getPeaksinIntensityOrder(Spectrum s, int number){
		s.setRanksOfPeaks();
		HashMap<Integer, Peak> peakMap = new HashMap<Integer, Peak>();
		for(Peak p:s){
			if(p.getRank() <= number)
				peakMap.put(p.getRank(), p);
		}
		
		return peakMap;
	}
	
	HashMap<Integer, Peak> getPeaksinIntensityOrderExcludingPrecursorIons(Spectrum s, int number){
		
		HashMap<Integer, Peak> peakMap = new HashMap<Integer, Peak>();
		
		for(Peak p:s){
			boolean isPrecursorIon = false;
			for(IonType ion : precursorIons){
				float toleranceInDa = tolerance.getToleranceAsDa(s.getParentMass());
				if(Math.abs(p.getMz() - ion.getMz(s.getParentMass())) <= toleranceInDa){
					isPrecursorIon = true;
					break;
				}
			}
			if(isPrecursorIon){
				p.setIntensity(0);
			}
		}
		
		s.setRanksOfPeaks();
		
		for(Peak p:s){
			if(p.getRank() <= number)
					peakMap.put(p.getRank(), p);
		}
		
	   //System.out.println(peakMap.size());
		return peakMap;
	}
	
	HashMap<Integer, Peak> getPeaksbyAlgorithm1(Spectrum s, int number){
		int n = Math.min(number, 10);
		int charge = s.getCharge();
		if(charge == 1) return null;
		
		ArrayList<Peak> sigPeaks = new ArrayList<Peak>();
		HashMap<Integer, Peak> peakMap = new HashMap<Integer, Peak>();
		HashMap<Peak, Peak> pp = new HashMap<Peak, Peak>();
		
		float alpha = 1, beta = 1;
		float maxIntensity = 0;
		
		Spectrum s2 = s.getCloneWithoutPeakList();
		
		s.setRanksOfPeaks();

		for(Peak p:s){
			if(p.getRank() <= n){
				sigPeaks.add(p);
				if(p.getRank() <= 1) maxIntensity = p.getIntensity();
			}
		}
		
		for(Peak p:sigPeaks){
			float score = 	0;
			
			if(p.getMz() > 0.55 * s.getParentMass()) score++;
			
			for(AminoAcid aa : aaSet){
				for(Peak p2 : s.getPeakListByMass(p.getMz() + aa.getMass()/(charge - 1), tolerance)){
					if(p2.getRank() <= n) score++;
				}
				
				for(Peak p2 : s.getPeakListByMass(p.getMz() - aa.getMass()/(charge - 1), tolerance)){
					if(p2.getRank() <= n) score++;
				}
				
				score *= alpha;
			}
			
			score -= getDerivativeNum(p, sigPeaks, s.getParentMass(), s.getCharge()) * beta;
			
			Peak newPeak = new Peak(p.getMz(), score + p.getIntensity() * Math.min(alpha, beta) * 0.9f / maxIntensity, p.getCharge());
			s2.add(newPeak);
			pp.put(newPeak, p);
			
		}
		
		s2.setRanksOfPeaks();
		
		
		for(Peak newPeak:s2){
			peakMap.put(newPeak.getRank(), pp.get(newPeak));
		}
		
		for(Peak p: s){
			if(p.getRank() > n && p.getRank() <= number)
				peakMap.put(p.getRank(), p);
		}
		
	//	System.out.println(getPeaksinIntensityOrder(s, number));
	//	System.out.println("*" + peakMap);
		return peakMap;

	}
	
	
	HashMap<Integer, Peak> getPeaksbyAlgorithm2(Spectrum s, int number){
		int n = number;
		int charge = s.getCharge();
		if(charge == 1) return null;
		
		ArrayList<Peak> sigPeaks = new ArrayList<Peak>();
		HashMap<Integer, Peak> peakMap = new HashMap<Integer, Peak>();
		HashMap<Peak, Peak> pp = new HashMap<Peak, Peak>();
		
		float alpha = 1, beta = 1;
		
		Spectrum s2 = s.getCloneWithoutPeakList();
		
		s.setRanksOfPeaks();

		for(Peak p:s){
			if(p.getRank() <= n){
				sigPeaks.add(p);
			}
		}
		
		for(Peak p:sigPeaks){
			float score = 	0;
			
			score += n - p.getRank();
			
			if(p.getMz() > 0.55 * s.getParentMass()) score++;
			
			for(AminoAcid aa : aaSet){
				for(Peak p2 : s.getPeakListByMass(p.getMz() + aa.getMass()/(charge - 1), tolerance)){
					if(p2.getRank() <= n) score++;
				}
				
				for(Peak p2 : s.getPeakListByMass(p.getMz() - aa.getMass()/(charge - 1), tolerance)){
					if(p2.getRank() <= n) score++;
				}
				
				score *= alpha;
			}
			
			score -= getDerivativeNum(p, sigPeaks, s.getParentMass(), s.getCharge()) * beta;
			
			Peak newPeak = new Peak(p.getMz(), score, p.getCharge());
			s2.add(newPeak);
			pp.put(newPeak, p);
			
		}
		
		s2.setRanksOfPeaks();
		
		
		for(Peak newPeak:s2){
			peakMap.put(newPeak.getRank(), pp.get(newPeak));
		}
		
		for(Peak p: s){
			if(p.getRank() > n && p.getRank() <= number)
				peakMap.put(p.getRank(), p);
		}
		
	//	System.out.println(getPeaksinIntensityOrder(s, number));
	//	System.out.println("*" + peakMap);
		return peakMap;

	}
	
	
	private int getDerivativeNum(Peak peak, ArrayList<Peak> sigPeaks, float precursorMass, int specCharge){
		int num = 0;
		
		for(IonType ion : precursorIons){
			float toleranceInDa = tolerance.getToleranceAsDa(precursorMass);
			if(Math.abs(peak.getMz() - ion.getMz(precursorMass)) <= toleranceInDa)
				num++;
		}

		for(Peak p : sigPeaks){
			float toleranceInDa = tolerance.getToleranceAsDa(p.getMz());
			float diff = p.getMz() - peak.getMz();
			
			if(peak.getRank() > p.getRank()){		
				// b and y
				if(Math.abs(precursorMass - peak.getMz() + specCharge * Composition.PROTON - p.getMz()) <= toleranceInDa)
					num++;
				// - h2o - nh3, + isotope
				if(Math.abs(diff - Composition.H2O) <= toleranceInDa)
					num++;
				if(Math.abs(diff - Composition.NH3) <= toleranceInDa)
					num++;
				if(Math.abs(-diff - Composition.ISOTOPE) <= toleranceInDa)
					num++;	
			}
			// y2 and y or b2 and b
			if(Math.abs(peak.getMz() * specCharge - Composition.PROTON - p.getMz() * (specCharge-1)) <= toleranceInDa)
				num++;
		}
		
		return num;
		
	}
	
	
}
