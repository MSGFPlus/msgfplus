package msgappednovo.parameters;

import java.util.ArrayList;

import msutil.Spectrum;

public class SpectrumParameter {
	private int specCharge;
	private int specMzRange;
	static private int specMzRangeNum = 5;
	
	private SpectrumParameter(int specCharge, int specMzRange){
		this.specCharge = specCharge;
		this.specMzRange = specMzRange;
	}
	
	public SpectrumParameter(Spectrum spec){
		this.specCharge = spec.getCharge();
		this.specMzRange = getSpecMzRange(spec);
	}
	
	public String toFileString() { return specCharge + " " + specMzRange;}
	
	public String toString(){
		return "Charge: " + specCharge + " Mz: " + specMzRange;
	}
	
	public boolean equals(Object o){
		if(this == o) return true;
    	if ( !(o instanceof SpectrumParameter) ) return false;
    	SpectrumParameter spar = (SpectrumParameter)o;
    	return this.specCharge == spar.specCharge && this.specMzRange == spar.specMzRange;
    
	}
	
	public int hashCode(){ return specCharge * (specMzRange + 1237);}
	
	public int getSpecCharge() { return specCharge; }
	public int getSpecMzRange() { return specMzRange; }
	
	static public int getSpecMzRange(Spectrum spec){ 
		int len = Math.round(spec.getParentMass() / 121.6f);
		//TODO divide using specMzRangeNum
		if(specMzRangeNum == 1) return 0;
		if(len <= 9) return 0;
		if(len <= 14) return 1;
		if(len <= 18) return 2;
		if(len <= 22) return 3;
		return 4;
	}
	
	static public void setSpecMzRangeNum(int n) {specMzRangeNum = n;}
	static public int getSpecMzRangeNum() {return specMzRangeNum;}
	
	static public SpectrumParameter parseFileString(String s){
		String[] token = s.split(" ");
		return new SpectrumParameter(Integer.parseInt(token[0]), Integer.parseInt(token[1]));
	}
	
	static public ArrayList<SpectrumParameter> getAllSpectrumParameters(int charge){
		ArrayList<SpectrumParameter> pars = new ArrayList<SpectrumParameter>();
	
		for(int r=0; r<5; r++){
			pars.add(new SpectrumParameter(charge, r));
		}
		
		return pars;
	}
}
