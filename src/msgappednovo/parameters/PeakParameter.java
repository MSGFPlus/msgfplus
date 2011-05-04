package msgappednovo.parameters;

import java.util.ArrayList;
import msutil.Peak;
import msutil.Spectrum;

public class PeakParameter {
	private static int peakIntensityRatioNum = 5 * 2;
	static private int groupNum = 10;
	static private int partitionNum = 4; // 0 => partition num varies
	static private Spectrum currentSpec = null;
	static private float currentHighestPeakIntensity = 0;
	
	private int basePeakGroup, basePeakPartition;
	
	private PeakParameter(int basePeakGroup, int basePeakPartition){
		this.basePeakGroup = basePeakGroup;
		this.basePeakPartition = basePeakPartition;
	}
	
	public PeakParameter(Peak bp, Spectrum spec, int iterationNum){
		if(currentSpec == null || currentSpec != spec){
			currentSpec = spec;
			currentHighestPeakIntensity = -1;
			for(Peak p : spec){
				if(p.getRank() == 1){
					currentHighestPeakIntensity = p.getIntensity();
					break;
				}
				//	currentHighestPeakIntensity = Math.max(currentHighestPeakIntensity, p.getIntensity());
			}
		}
		
		if(groupNum == 1) basePeakGroup = 0;
		else if(iterationNum == 0)
			basePeakGroup = Math.min(groupNum - 1, bp.getRank()/(100/(groupNum-1)));
		else
			basePeakGroup = groupNum - 1 - Math.round((groupNum-2) * bp.getIntensity()/currentHighestPeakIntensity);

		if(partitionNum <= 0){
			basePeakPartition = 1;
			for(int charge = spec.getCharge(); charge >=2; charge--){
				if(bp.getMz() < maxMzWith(charge, spec)){
					basePeakPartition = charge;
					break;
				}
			}
		}else{
			float partitionSize = spec.getParentMass() / partitionNum;
			basePeakPartition =  (int)Math.min((bp.getMz() / partitionSize), partitionNum-1);
		}
			//if(bp.getMz() > spec.getParentMass()/4*3) basePeakPartition = spec.getCharge() + 1;
			//float partitionSize = spec.getParentMass() / partitionNum;
			//basePeakPartition =  (int)Math.min((bp.getMz() / partitionSize) + 1, partitionNum);
		//}
	}
	
	public String toFileString(){
		return basePeakGroup + " " + basePeakPartition;
	}
	
	public String toString(){
		return "Group: " + basePeakGroup + " Partition: " + basePeakPartition;
	}
	
	public int getBasePeakGroupNum() { return basePeakGroup; }
	public int getBasePeakPartitionNum() { return basePeakPartition; }
	
	public boolean equals(Object o){
		if(this == o) return true;
    	if ( !(o instanceof PeakParameter) ) return false;
    	PeakParameter ppar = (PeakParameter)o;
    	
    	return this.basePeakGroup == ppar.basePeakGroup 
    		&& this.basePeakPartition == ppar.basePeakPartition;
	}
	
	public int hashCode(){
		return basePeakGroup * (basePeakPartition + 12373);
	}
	
	static public void setGroupNum(int g) { groupNum = g; }
	
	static public void setPartitionNum(int p) {partitionNum = p;}
	
	static public int getMaxGroupNum() {return groupNum;}
	
	static public PeakParameter parseFileString(String s){
		String[] token = s.split(" ");
		return new PeakParameter(Integer.parseInt(token[0]), Integer.parseInt(token[1]));
	}
	
	static public float maxMzWith(int charge, Spectrum spec){		
		if(charge == 1) return Float.MAX_VALUE;
		
		return (float) (spec.getParentMass()/charge);
	}
	
	static public ArrayList<PeakParameter> getAllBasePeakParameters(int charge){
		ArrayList<PeakParameter> pars = new ArrayList<PeakParameter>();
		float maxPartitionNum = partitionNum-1;
		if(partitionNum <= 0) maxPartitionNum = charge;
		
		for(int p=0; p<=maxPartitionNum; p++){
			for(int g=0; g<groupNum; g++){
				pars.add(new PeakParameter(g, p));
			}
		}
		return pars;
	}
	
	static public void setPeakIntensityRatioNum(int n){ peakIntensityRatioNum = n;}
	
	static public int getPeakIntensityRatioNum(Peak bp, Peak cp, Spectrum spec){
		float bpi = bp.getIntensity();
		float cpi = cp.getIntensity();
		
		if(bpi == 0) bpi = Float.MIN_VALUE;
		if(cpi == 0) cpi = Float.MIN_VALUE;
		
		if(peakIntensityRatioNum == 1) return 0; 
		
		if(spec.getPrecursorPeak().equals(cp)) return -peakIntensityRatioNum/2;
		
		float x = bpi/cpi;
		
		if(x>=1){
			return  -(int) Math.floor((1-1/x)*(peakIntensityRatioNum/2));
			}
		else{
			return  (int) Math.floor((1-x)*(peakIntensityRatioNum/2)) + 1;
		}
	}
	
	static public void refresh() {currentSpec = null;}
	
	static public ArrayList<Integer> getAllPeakIntensityRatioNums(){
		ArrayList<Integer> ratios = new ArrayList<Integer>();

		for(int r = -peakIntensityRatioNum/2; r<= peakIntensityRatioNum/2; r++)
			ratios.add(r);		
		return ratios;
	}

	
}
