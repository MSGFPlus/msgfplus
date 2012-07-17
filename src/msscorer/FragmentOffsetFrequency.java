package msscorer;

import msutil.IonType;

public class FragmentOffsetFrequency implements Comparable<FragmentOffsetFrequency> {
	public FragmentOffsetFrequency(IonType ionType, float frequency) {
		super();
		this.ionType = ionType;
		this.frequency = frequency;
	}
	public IonType getIonType() {
		return ionType;
	}
	public void setIonType(IonType ionType) {
		this.ionType = ionType;
	}
	public float getFrequency() {
		return frequency;
	}
	public void setFrequency(float probability) {
		this.frequency = probability;
	}
	public int compareTo(FragmentOffsetFrequency o) {
		if(this.frequency > o.frequency)
			return 1;
		else if(this.frequency == o.frequency)
			return 0;
		else
			return -1;
	}
	
	private IonType ionType;
	private float frequency;
}
