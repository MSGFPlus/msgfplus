package msutil;

import java.util.ArrayList;
import java.util.Iterator;

public class SpecKey extends Pair<Integer, Integer> {

	public SpecKey(int scanNum, int charge) {
		super(scanNum, charge);
	}
	
	public int getScanNum()
	{
		return super.getFirst();
	}
	
	public int getCharge()
	{
		return super.getSecond();
	}
	
	public String getSpecKeyString()
	{
		return getScanNum()+":"+getCharge();
	}

	public static SpecKey getSpecKey(String specKeyString)
	{
		String[] token = specKeyString.split(":");
		return new SpecKey(Integer.parseInt(token[0]), Integer.parseInt(token[1]));
	}
	
	public static ArrayList<SpecKey> getSpecKeyList(Iterator<Spectrum> itr, int minCharge, int maxCharge)
	{
		ArrayList<SpecKey> specKeyList = new ArrayList<SpecKey>();
		
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			int scanNum = spec.getScanNum();
			int charge = spec.getCharge();
			if(charge == 0)
			{
				for(int c=minCharge; c<=maxCharge; c++)
					specKeyList.add(new SpecKey(scanNum, c));
			}
			else if(charge > 0)
			{
				specKeyList.add(new SpecKey(scanNum, charge));
			}
		}
		return specKeyList;
	}
}
