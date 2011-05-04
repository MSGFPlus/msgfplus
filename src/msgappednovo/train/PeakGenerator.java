package msgappednovo.train;

import java.util.ArrayList;
import java.util.HashMap;

import msgf.Tolerance;
import msutil.Composition;
import msutil.IonType;
import msutil.Peak;
import msutil.Peptide;
import msutil.Spectrum;

public class PeakGenerator {
	static final float VirtualPeakIntensity = Float.POSITIVE_INFINITY;
	private Peptide pep;
	private HashMap<IonType, ArrayList<Peak>> ionTheoreticalBasePeakMap = null; // theoretical base peaks
	
	public PeakGenerator(Spectrum spec){
		pep = spec.getAnnotation();
	}
	
	public boolean isExplainedBy(Peak p, IonType ion, Tolerance tol, Tolerance pmtol){
		boolean isPrecursorIon = (ion instanceof IonType.PrecursorIon);
		Tolerance tmptol = isPrecursorIon? pmtol : tol;
		for(Peak peak : getTheoreticalBasePeaks(ion)){
			float diff = Math.abs(p.getMz() - peak.getMz());
			if(diff <= tmptol.getToleranceAsDa(p.getMass())/ (isPrecursorIon ? 1 :ion.getCharge())){
			//	if(ion.getCharge() == 2){
			//		System.out.println(ion + "\t" + peak.getMz() + "\t" + p.getMz()+"\t" + this.pep);
			//	}
				return true;
			}
		}
		return false;
	}
	
	ArrayList<Peak> getTheoreticalPrefixBasePeaks(int charge){ return getTheoreticalPrefixBasePeaks(charge, 0); }
	
	ArrayList<Peak> getTheoreticalPrefixBasePeaks(int charge, float massOffset){
		ArrayList<Peak> peaks = new ArrayList<Peak>();
		
		if(pep == null) return null;
		
		for(float m : pep.getPRMMasses(true, 0)){
			peaks.add(new Peak((m + massOffset)/charge, VirtualPeakIntensity, charge));
		}
		
		return peaks;
		
	}
	
	ArrayList<Peak> getTheoreticalSuffixBasePeaks(int charge){ return getTheoreticalSuffixBasePeaks(charge, 0); }
	
	ArrayList<Peak> getTheoreticalSuffixBasePeaks(int charge, float massOffset){
		ArrayList<Peak> peaks = new ArrayList<Peak>();
		
		if(pep == null) return null;
		
		for(float m : pep.getPRMMasses(false, 0)){
			peaks.add(new Peak((m + massOffset)/charge, VirtualPeakIntensity, charge));
		}
		
		return peaks;
		
	}
	
	ArrayList<Peak> getTheoreticalInternalBasePeaks(int charge){ return getTheoreticalInternalBasePeaks(charge, 0); }
	
	ArrayList<Peak> getTheoreticalInternalBasePeaks(int charge, float massOffset){
		ArrayList<Peak> peaks = new ArrayList<Peak>();
		
		if(pep == null) return null;
		
		for(int i=1; i<pep.size()-1; i++){
			for(int j=i+1; j<pep.size()-1;j++){
				float m = pep.getMass(i, j);
				peaks.add(new Peak((m + massOffset)/charge, VirtualPeakIntensity, charge));
			}
		}

		return peaks;
		
	}
	
	Peak getTheoreticalPrecursorBasePeak(int charge){ return getTheoreticalPrecursorBasePeak(charge, 0); }
	
	Peak getTheoreticalPrecursorBasePeak(int charge, float massOffset) { 
		if(pep == null) return null;
		return new Peak((float) ((pep.getParentMass() + massOffset)/charge + Composition.PROTON), VirtualPeakIntensity, charge); 
	}

	ArrayList<Peak> getTheoreticalBasePeaks(IonType ion){
		ArrayList<Peak> peaks;
		
		if(pep == null) return null;
		
		if(ionTheoreticalBasePeakMap == null)
			ionTheoreticalBasePeakMap = new HashMap<IonType, ArrayList<Peak>>();
		
		if((peaks = ionTheoreticalBasePeakMap.get(ion)) != null)
			return peaks;
		else peaks = new ArrayList<Peak>();
		
		if(ion instanceof IonType.PrefixIon){
			peaks = this.getTheoreticalPrefixBasePeaks(ion.getCharge(), ion.getOffset() * ion.getCharge());
		}else if(ion instanceof IonType.SuffixIon){
			peaks = this.getTheoreticalSuffixBasePeaks(ion.getCharge(), ion.getOffset() * ion.getCharge());
		}else if(ion instanceof IonType.PrecursorIon){
			peaks.add(this.getTheoreticalPrecursorBasePeak(ion.getCharge(), ion.getOffset() * ion.getCharge()));
		}else if(ion instanceof IonType.InternalIon){
			peaks = this.getTheoreticalInternalBasePeaks(ion.getCharge(), ion.getOffset() * ion.getCharge());
		}
		
		ionTheoreticalBasePeakMap.put(ion, peaks);
		
		return peaks;
	}
	
	Peak getTheoreticalComplementaryBasePeak(Peak p, int charge){
		if(pep == null) return null;
		return new Peak((float) (pep.getParentMass()/charge - p.getMz() + 2 * Composition.PROTON), VirtualPeakIntensity, charge);
	}
	
	static public Peak getChargeChangedBasePeak(Peak p, int chargeBefore, int chargeOffset){
		float newmz;
		int chargeAfter = chargeOffset + chargeBefore;
		
		if(chargeBefore == chargeAfter){
			newmz = p.getMz();
		}else{
			float mass = (float) ((p.getMz() - Composition.PROTON) * chargeBefore);
			newmz = (float) (mass/chargeAfter + Composition.PROTON);
		}
		
		return new Peak(newmz, VirtualPeakIntensity, chargeAfter);
	}
	
	static public Peak getComplementaryBasePeak(Peak p, int charge, Spectrum spec){
		return new Peak((float) (spec.getParentMass()/charge - p.getMz() + 2 * Composition.PROTON), VirtualPeakIntensity, charge);
	}
	
	static public Peak getComplementaryBasePeak(Peak p, int charge, float pm){
		return new Peak((float) (pm/charge - p.getMz() + 2 * Composition.PROTON), VirtualPeakIntensity, charge);
	}
	
	static public float getPrefixMass(Peak p, IonType ion, Spectrum spec){
		float mass = ion.getMass(p.getMz());
		
		if(ion instanceof IonType.SuffixIon){
			mass = (float) (spec.getParentMass() - Composition.H2O - mass);
		}
		
		return mass;
	}
}
