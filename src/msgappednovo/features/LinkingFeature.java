package msgappednovo.features;

import java.util.ArrayList;
import java.util.HashMap;

import msgappednovo.parameters.PeakParameter;
import msgappednovo.parameters.SpectrumParameter;
import msgf.Tolerance;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Peak;
import msutil.Spectrum;

public class LinkingFeature extends Feature{
	private int peakIntensityRatio;
	private AminoAcidSet aaSet;
	
	static private HashMap<Peak, HashMap<Integer, ArrayList<Integer>>> peakRatioMap = null;
	static private Spectrum currentSpec = null;
	static private int currentIterationNum = -1;
	
	public LinkingFeature(SpectrumParameter spar,
			PeakParameter ppar, int basePeakCharge, boolean isPresent,
			int iterationNum, int peakIntensityRatio, AminoAcidSet aaSet) {
		super(spar, ppar, basePeakCharge, isPresent, iterationNum);
		this.peakIntensityRatio = peakIntensityRatio;
		this.aaSet = aaSet;
	}

	private Peak currentPeak = null;
	private boolean currentHold = true;
	@Override
	public boolean holdsFor(Peak bp, Spectrum spec, Tolerance tol, Tolerance pmtol, int iterationNum) {
		if(spec.equals(currentSpec)&& bp == currentPeak){
			return currentHold;
		}
		
		boolean hold = false;

		if(super.holdsFor(bp, spec, iterationNum)){	
			if(!spec.equals(currentSpec) || currentIterationNum != iterationNum){
				currentSpec = spec; currentIterationNum = iterationNum;
				peakRatioMap = new HashMap<Peak, HashMap<Integer, ArrayList<Integer>>>();
			}
			
			if(!peakRatioMap.containsKey(bp))
				peakRatioMap.put(bp, new HashMap<Integer, ArrayList<Integer>>());
			HashMap<Integer, ArrayList<Integer>> r = peakRatioMap.get(bp);
				
			if(!r.containsKey(bp.getCharge()))
				r.put(bp.getCharge(),  new ArrayList<Integer>());
			
			ArrayList<Integer> ratios = r.get(bp.getCharge());
			
			ratios = new ArrayList<Integer>();
			for(AminoAcid aa : aaSet){
				for(int i=0; i<2; i++){
					float mz = bp.getMz() + aa.getMass()/bp.getCharge() * (i == 0 ? 1 : -1);
					float t = tol.getToleranceAsDa(mz)/this.getBasePeakCharge();
					for(Peak cp : spec.getPeakListByMassRange(mz - t, mz + t)){
						ratios.add(PeakParameter.getPeakIntensityRatioNum(bp, cp, spec));
					}
				}
			}
		
			
			hold = ratios.contains(this.peakIntensityRatio) == this.isPresent();	
		}
		currentPeak = bp;
		currentHold = hold;
		return hold;
	}

	@Override
	public boolean equals(Object o) {
		if(this == o) return true;
    	if ( !(o instanceof LinkingFeature) ) return false;
    	LinkingFeature con = (LinkingFeature)o;
    	
    	return this.peakIntensityRatio == con.peakIntensityRatio && super.equals(con);
	}


	@Override
	public int hashCode() {
		return this.getBasePeakParameter().hashCode() * this.getSpectrumParameter().hashCode() * this.getBasePeakCharge() 
		* (this.getIterationNum() + 531) * (this.peakIntensityRatio + 12773) * new Boolean(this.isPresent()).hashCode();
	}

	@Override
	public String toString() {
		return "Bridging - Iteration: " + this.getIterationNum() + " " + this.getSpectrumParameter() + " " + this.getBasePeakParameter() + " PeakCharge: " + this.getBasePeakCharge() + " Ratio: " + this.peakIntensityRatio + " Present: " + this.isPresent();
	}

	@Override
	public String toFileString() {
		return "B\t" + this.getIterationNum() + "\t" + this.getSpectrumParameter().toFileString() + "\t" + this.getBasePeakParameter().toFileString()+"\t" +this.getBasePeakCharge() + "\t" + this.peakIntensityRatio + "\t" + (this.isPresent()? "p" : "a");
	}
	
	static public LinkingFeature parseFileString(String s, AminoAcidSet aaSet){
		String[] token = s.split("\t");
		int iterationNum = Integer.parseInt(token[1]);
		SpectrumParameter spar = SpectrumParameter.parseFileString(token[2]);
		PeakParameter ppar = PeakParameter.parseFileString(token[3]);
		int basePeakCharge = Integer.parseInt(token[4]);
		int peakIntensityRatio = Integer.parseInt(token[5]);
		boolean isPresent = token[6].equals("p");
		
		return new LinkingFeature(spar, ppar, basePeakCharge, isPresent, iterationNum, peakIntensityRatio, aaSet);
	}

}
