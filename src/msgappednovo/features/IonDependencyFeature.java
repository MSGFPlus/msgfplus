package msgappednovo.features;

import msgappednovo.parameters.PeakParameter;
import msgappednovo.parameters.SpectrumParameter;
import msgappednovo.train.InterPeakOffset;
import msgf.Tolerance;
import msutil.Peak;
import msutil.Spectrum;

public class IonDependencyFeature extends Feature{
	private int peakIntensityRatio;
	private InterPeakOffset gof;
	
	public IonDependencyFeature(SpectrumParameter spar,
			PeakParameter ppar, int basePeakCharge, boolean isPresent,
			int iterationNum, int peakIntensityRatio, InterPeakOffset gof) {
		super(spar, ppar, basePeakCharge, isPresent, iterationNum);
		this.peakIntensityRatio = peakIntensityRatio;
		assert(basePeakCharge == gof.getBaseCharge());
		this.gof = gof;
	}


	private Peak currentPeak = null;
	private boolean currentHold = true;
	private Spectrum currentSpec = null;
	
	@Override
	public boolean holdsFor(Peak bp, Spectrum spec, Tolerance tol, Tolerance pmtol, int iterationNum) {
		if(spec.equals(currentSpec) && bp == currentPeak){
			return currentHold;
		}
		
		boolean hold = false;

		if(super.holdsFor(bp, spec, iterationNum)){
			boolean match = false;
		
			for(Peak cp : gof.getMatchingPeaks(bp, spec, tol, pmtol)){
				if(peakIntensityRatio == PeakParameter.getPeakIntensityRatioNum(bp, cp, spec)){
					match = true;
					break;
				}
			}
			hold = match == this.isPresent();	
		}
		
		currentSpec = spec;
		currentHold = hold;
		currentPeak = bp;
		
		return hold;
	}

	@Override
	public boolean equals(Object o) {
		if(this == o) return true;
    	if ( !(o instanceof IonDependencyFeature) ) return false;
    	IonDependencyFeature con = (IonDependencyFeature)o;
    	
    	return this.peakIntensityRatio == con.peakIntensityRatio && super.equals(con) && this.gof.equals(con.gof);
	}

	@Override
	public int hashCode() {
		return this.getBasePeakParameter().hashCode() * this.getSpectrumParameter().hashCode() * this.getBasePeakCharge() 
		* (this.getIterationNum() + 7681) * (this.peakIntensityRatio + 48701)  * new Boolean(this.isPresent()).hashCode()
		* this.gof.hashCode();
	}

	@Override
	public String toString() {
		return "GOF - Iteration: " + this.getIterationNum() + " " + this.getSpectrumParameter() + " " + this.getBasePeakParameter() + " PeakCharge: " + this.getBasePeakCharge() + " Ratio: " + this.peakIntensityRatio + " Present: " + this.isPresent() + " GOF: " + this.gof;
	}
	
	@Override
	public String toFileString() {
		return "G\t" + this.getIterationNum() + "\t" + this.getSpectrumParameter().toFileString() + "\t" + this.getBasePeakParameter().toFileString() +"\t" +this.getBasePeakCharge() + "\t" + this.peakIntensityRatio + "\t" + (this.isPresent()? "p" : "a") + "\t" + this.gof.toFileString();
	}
	
	public boolean isComplementary() {return gof.isComplementary();}
	public int getChargeOffset() {return gof.getChargeOffset();}
	
	static public IonDependencyFeature parseFileString(String s){
		String[] token = s.split("\t");
		int iterationNum = Integer.parseInt(token[1]);
		SpectrumParameter spar = SpectrumParameter.parseFileString(token[2]);
		PeakParameter ppar = PeakParameter.parseFileString(token[3]);
		int basePeakCharge = Integer.parseInt(token[4]);
		int peakIntensityRatio = Integer.parseInt(token[5]);
		boolean isPresent = token[6].equals("p");
		InterPeakOffset gof = InterPeakOffset.parseFileString(token[7]);
		
		return new IonDependencyFeature(spar, ppar, basePeakCharge, isPresent, iterationNum, peakIntensityRatio, gof);
	}
	
}
