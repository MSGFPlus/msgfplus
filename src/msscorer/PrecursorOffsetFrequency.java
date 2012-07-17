package msscorer;

import java.util.ArrayList;

import msgf.Tolerance;

public class PrecursorOffsetFrequency implements Comparable<PrecursorOffsetFrequency> {
	public PrecursorOffsetFrequency(int reducedCharge, float offset, float frequency) {
		super();
		this.reducedCharge = reducedCharge;
		this.offset = offset;
		this.frequency = frequency;
		this.tolerance = new Tolerance(0.5f);
	}
	
	public PrecursorOffsetFrequency tolerance(Tolerance tolerance)
	{
		this.tolerance = tolerance;
		return this;
	}
	
	public int getReducedCharge() {
		return reducedCharge;
	}
	public void setReducedCharge(int reducedCharge) {
		this.reducedCharge = reducedCharge;
	}
	public float getOffset() {
		return offset;
	}
	public void setOffset(float offset) {
		this.offset = offset;
	}
	public float getFrequency() {
		return frequency;
	}
	public void setFrequency(float probability) {
		this.frequency = probability;
	}
	public Tolerance getTolerance()
	{
		return tolerance;
	}
	
	private int reducedCharge;
	private float offset;
	private float frequency;
	private Tolerance tolerance;
	
	public int compareTo(PrecursorOffsetFrequency o) {
		return new Float(this.frequency).compareTo(new Float(o.frequency));
	}
	
	public static ArrayList<PrecursorOffsetFrequency> getClusteredOFF(ArrayList<PrecursorOffsetFrequency> offList, float granularity)
	{
		ArrayList<PrecursorOffsetFrequency> clusteredOFF = new ArrayList<PrecursorOffsetFrequency>();
		if(offList == null)
			return null;
		else if(offList.size() == 0)
			return clusteredOFF;
		
		PrecursorOffsetFrequency prevOFF = offList.get(0);
		int clusterStartIndex = 0;
		float clusterFreq = prevOFF.getFrequency();
		int reducedCharge = prevOFF.getReducedCharge();
		
		for(int i=1; i<offList.size(); i++)
		{
			PrecursorOffsetFrequency off = offList.get(i);
			if(Math.abs(off.getOffset()-prevOFF.getOffset()-granularity) < granularity*0.1f)
			{
				clusterFreq += off.getFrequency();
			}
			else
			{
				float offset = (offList.get(clusterStartIndex).getOffset()+offList.get(i-1).getOffset())/2;
				float tolDa = granularity/2*(i-clusterStartIndex);
				clusteredOFF.add(new PrecursorOffsetFrequency(reducedCharge, offset, clusterFreq).tolerance(new Tolerance(tolDa)));
				clusterStartIndex = i;
				clusterFreq = off.getFrequency();
			}
			prevOFF = off;
		}
		float offset = offList.get(clusterStartIndex).getOffset()+offList.get(offList.size()-1).getOffset()/2;
		float tolDa = granularity/2*(offList.size()-clusterStartIndex);
		clusteredOFF.add(new PrecursorOffsetFrequency(reducedCharge, offset, clusterFreq).tolerance(new Tolerance(tolDa)));
		
//		for(PrecursorOffsetFrequency off : clusteredOFF)
//			System.out.println(off.getReducedCharge()+"\t"+off.getOffset()+"\t"+off.getFrequency()+"\t"+off.getTolerance().toString());
		return clusteredOFF;
	}
}
