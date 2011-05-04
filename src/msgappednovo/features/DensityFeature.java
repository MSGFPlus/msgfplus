package msgappednovo.features;

import java.util.ArrayList;
import java.util.HashMap;

import msgappednovo.parameters.PeakParameter;
import msgappednovo.parameters.SpectrumParameter;
import msgf.Tolerance;
import msutil.AminoAcid;
import msutil.Composition;
import msutil.Peak;
import msutil.Spectrum;

public class DensityFeature extends Feature{
	private int peakIntensityRatio;
	
	// variables for speed up
	static private HashMap<Peak, HashMap<Integer, ArrayList<Integer>>> peakRatioMap = null;
	static private Spectrum currentSpec = null;
	static private int currentIterationNum = -1;
	
	public DensityFeature(SpectrumParameter spar,
			PeakParameter ppar, int basePeakCharge, boolean isPresent,
			int iterationNum, int peakIntensityRatio) {
		super(spar, ppar, basePeakCharge, isPresent, iterationNum);
		this.peakIntensityRatio = peakIntensityRatio;
	}

	private Peak currentPeak = null;
	private boolean currentHold = true;
	@Override
	public boolean holdsFor(Peak bp, Spectrum spec, Tolerance tol, Tolerance pmtol, int iterationNum) {
		if(spec.equals(currentSpec)&& bp ==  currentPeak){
			return currentHold;
		}
		boolean hold = false;
		
		if(super.holdsFor(bp, spec, iterationNum)){
			if(!spec.equals(currentSpec) || currentIterationNum != iterationNum){
				currentSpec = spec; currentIterationNum = iterationNum;
				peakRatioMap = new HashMap<Peak, HashMap<Integer, ArrayList<Integer>>>();
			}
			
			if(!peakRatioMap.containsKey(bp)) peakRatioMap.put(bp, new HashMap<Integer, ArrayList<Integer>>());
			HashMap<Integer, ArrayList<Integer>> r = peakRatioMap.get(bp);
			
			if(!r.containsKey(bp.getCharge())) r.put(bp.getCharge(), new ArrayList<Integer>());
			ArrayList<Integer> ratios = r.get(bp.getCharge());
			
			if(ratios.isEmpty()){
				float mz = bp.getMz();
				float min = (float) (Composition.ISOTOPE + tol.getToleranceAsDa(mz)/this.getBasePeakCharge());
				float max = AminoAcid.getStandardAminoAcid('G').getMass()- tol.getToleranceAsDa(mz)/this.getBasePeakCharge();
				
				for(int i=0; i<2; i++){	
					ArrayList<Peak> peakList;
					
					if(i == 0) peakList = spec.getPeakListByMassRange(mz + min, mz + max);
					else peakList = spec.getPeakListByMassRange(mz - max, mz - min);
					
					for(Peak cp : peakList){
						ratios.add(PeakParameter.getPeakIntensityRatioNum(bp, cp, spec));
					}
				}
			}
			
			hold =  ratios.contains(this.peakIntensityRatio) == this.isPresent();	
		}
		return hold;
	}

	@Override
	public boolean equals(Object o) {
		if(this == o) return true;
    	if ( !(o instanceof DensityFeature) ) return false;
    	DensityFeature con = (DensityFeature)o;
    	
    	return this.peakIntensityRatio == con.peakIntensityRatio && super.equals(con);
	}

	@Override
	public int hashCode() {
		return -this.getBasePeakParameter().hashCode() * this.getSpectrumParameter().hashCode() * this.getBasePeakCharge() 
		* (this.getIterationNum() + 9641) * (this.peakIntensityRatio + 16543)  * new Boolean(this.isPresent()).hashCode();
	}

	@Override
	public String toString() {
		return "Density - Iteration: " + this.getIterationNum() + " " + this.getSpectrumParameter() + " " + this.getBasePeakParameter() + " PeakCharge: " + this.getBasePeakCharge() + " Ratio: " + this.peakIntensityRatio + " Present: " + this.isPresent();
	}

	@Override
	public String toFileString() {
		return "D\t" + this.getIterationNum() + "\t" + this.getSpectrumParameter().toFileString() + "\t" + this.getBasePeakParameter().toFileString() +"\t" +this.getBasePeakCharge() + "\t" + this.peakIntensityRatio + "\t" + (this.isPresent()? "p" : "a");
	}

	static public DensityFeature parseFileString(String s){
		String[] token = s.split("\t");
		int iterationNum = Integer.parseInt(token[1]);
		SpectrumParameter spar = SpectrumParameter.parseFileString(token[2]);
		PeakParameter ppar = PeakParameter.parseFileString(token[3]);
		int basePeakCharge = Integer.parseInt(token[4]);
		int peakIntensityRatio = Integer.parseInt(token[5]);
		boolean isPresent = token[6].equals("p");
		
		return new DensityFeature(spar, ppar, basePeakCharge, isPresent, iterationNum, peakIntensityRatio);
	}
	
}
