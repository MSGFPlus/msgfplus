package edu.ucsd.msjava.msscorer;

public class Partition implements Comparable<Partition> {
	private int charge;
	private float parentMass;
	private int segIndex;

	public Partition(int charge, float parentMass, int segIndex) {
		super();
		this.charge = charge;
		this.parentMass = parentMass;
		this.segIndex = segIndex;
	}
	public int getCharge() {
		return charge;
	}

	public void setCharge(int charge) {
		this.charge = charge;
	}

	public float getParentMass() {
		return parentMass;
	}

	public void setParentMass(float parentMass) {
		this.parentMass = parentMass;
	}
	public int getSegNum() {
		return segIndex;
	}
	public void setPosIndex(int posIndex) {
		this.segIndex = posIndex;
	}
	public int compareTo(Partition o) 
	{
		if(charge < o.charge)
			return -1;
		else if(charge > o.charge)
			return 1;
		else 
		{
			if(segIndex < o.segIndex)
				return -1;
			else if(segIndex > o.segIndex)
				return 1;
			else
			{
				if(parentMass < o.parentMass)
					return -1;
				else if(parentMass > o.parentMass)
					return 1;
				else
					return 0;
			}
		}
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if(obj instanceof Partition)
		{
			Partition o = (Partition)obj;
			if(charge == o.charge && parentMass == o.parentMass && segIndex == o.segIndex)
				return true;
		}
		return false;
	}
	
	@Override
	public int hashCode()
	{
		return new Float(parentMass).hashCode()+charge*10+segIndex;
	}
}	
