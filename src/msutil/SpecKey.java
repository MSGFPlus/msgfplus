package msutil;

import java.util.ArrayList;
import java.util.Iterator;

import parser.MzXMLSpectraIterator;

public class SpecKey extends Pair<Integer, Integer> {

	private ArrayList<Integer> scanNumList;

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
	
	public void addScanNum(int scanNum)
	{
		if(scanNumList == null)
		{
			scanNumList = new ArrayList<Integer>();
			scanNumList.add(super.getFirst());
		}
		scanNumList.add(scanNum);
	}
	
	public ArrayList<Integer> getScanNumList()
	{
		return scanNumList;
	}
	
	public static ArrayList<SpecKey> getSpecKeyList(Iterator<Spectrum> itr, int minCharge, int maxCharge, ActivationMethod activationMethod)
	{
		if(activationMethod == ActivationMethod.FUSION)
			return getFusedSpecKeyList(itr, minCharge, maxCharge);
		
		ArrayList<SpecKey> specKeyList = new ArrayList<SpecKey>();
		
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			int scanNum = spec.getScanNum();
			int charge = spec.getCharge();
			
			if(activationMethod != null && spec.getActivationMethod() != null && spec.getActivationMethod() != activationMethod)
				continue;

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
	
	public static ArrayList<SpecKey> getFusedSpecKeyList(Iterator<Spectrum> itr, int minCharge, int maxCharge)
	{
		ArrayList<SpecKey> specKeyList = new ArrayList<SpecKey>();
		
		int prevScanNum = Integer.MIN_VALUE;
		int prevCharge = Integer.MIN_VALUE;
		float previousPrecursorMz = Float.MIN_VALUE;
		
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			int scanNum = spec.getScanNum();
			int charge = spec.getCharge();
			float precursorMz = spec.getPrecursorPeak().getMz();
			if(spec.getActivationMethod() == null)
			{
				System.out.println("Error: activation method is not available: Scan=" + spec.getScanNum()+", PrecursorMz=" + precursorMz);
				System.exit(-1);
			}
			
			if(scanNum == prevScanNum+1 && charge == prevCharge && precursorMz == previousPrecursorMz)
			{
				if(charge == 0)
				{
					for(int i=0; i<maxCharge-minCharge; i++)
						specKeyList.get(specKeyList.size()-1-i).addScanNum(scanNum);
				}
				else if(charge > 0)
				{
					specKeyList.get(specKeyList.size()-1).addScanNum(scanNum);
				}
			}
			else
			{
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
			prevScanNum = scanNum;
			prevCharge = charge;
			previousPrecursorMz = precursorMz;
		}
		return specKeyList;
	}
	
	public static void main(String argv[]) throws Exception
	{
		test();
	}
	
	public static void test() throws Exception
	{
		String fileName = "/home/sangtaekim/Research/Data/HeckRevision/CIDETDPairs/mzXML/090121_NM_Trypsin_20.mzXML";
		int minCharge = 2, maxCharge = 3;
		MzXMLSpectraIterator itr = new MzXMLSpectraIterator(fileName);
		ArrayList<SpecKey> list = SpecKey.getFusedSpecKeyList(itr, minCharge, maxCharge);
		for(SpecKey specKey : list)
		{
			if(specKey.getScanNumList() == null || specKey.getScanNumList().size() != 2)
				System.out.println(specKey.getSpecKeyString()+"\t"+specKey.getScanNumList());
		}
		System.out.println("Size: " + list.size());
	}
	
}
