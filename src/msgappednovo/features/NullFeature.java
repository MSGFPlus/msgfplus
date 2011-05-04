package msgappednovo.features;


import msgappednovo.parameters.PeakParameter;
import msgappednovo.parameters.SpectrumParameter;
import msgf.Tolerance;
import msutil.Peak;
import msutil.Spectrum;

public class NullFeature extends Feature{

	public NullFeature(SpectrumParameter spar,
			PeakParameter ppar, int basePeakCharge,
			int iterationnum) {
		super(spar, ppar, basePeakCharge, true, iterationnum);
	}

	@Override
	public boolean holdsFor(Peak bp, Spectrum spec, Tolerance tol,  Tolerance pmtol, int iterationNum) {
		return super.holdsFor(bp, spec, iterationNum);
	}

	public boolean equals(Object o) {
		if(this == o) return true;
    	if ( !(o instanceof NullFeature) ) return false;
    	NullFeature con = (NullFeature)o;
    	
    	return super.equals(con);
	}

	@Override
	public int hashCode() {
		return this.getBasePeakParameter().hashCode() * this.getSpectrumParameter().hashCode() * this.getBasePeakCharge() 
			* (this.getIterationNum() + 2451);
	}

	@Override
	public String toString() {
		return "Null - Iteration: " + this.getIterationNum() + " " + this.getSpectrumParameter() + " " + this.getBasePeakParameter() + " PeakCharge: " + this.getBasePeakCharge();
	}

	@Override
	public String toFileString() {
		return "N\t" + this.getIterationNum() + "\t" + this.getSpectrumParameter().toFileString() + "\t" + this.getBasePeakParameter().toFileString() +"\t" +this.getBasePeakCharge();
	}
	
	static public NullFeature parseFileString(String s){
		String[] token = s.split("\t");
		int iterationNum = Integer.parseInt(token[1]);
		SpectrumParameter spar = SpectrumParameter.parseFileString(token[2]);
		PeakParameter ppar = PeakParameter.parseFileString(token[3]);
		int basePeakCharge = Integer.parseInt(token[4]);
		
		return new NullFeature(spar,	ppar, basePeakCharge, iterationNum);
	}

}
