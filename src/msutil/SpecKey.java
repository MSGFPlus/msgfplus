package msutil;

import java.util.ArrayList;
import java.util.Iterator;

import parser.MzXMLSpectraIterator;
import sequences.Constants;

public class SpecKey extends Pair<Integer, Integer> {

	private ArrayList<Integer> specIndexList;

	public SpecKey(int specIndex, int charge) {
		super(specIndex, charge);
	}
	
	public int getSpecIndex()
	{
		return super.getFirst();
	}
	
	public int getCharge()
	{
		return super.getSecond();
	}
	
	public String getSpecKeyString()
	{
		return getSpecIndex()+":"+getCharge();
	}

	public static SpecKey getSpecKey(String specKeyString)
	{
		String[] token = specKeyString.split(":");
		return new SpecKey(Integer.parseInt(token[0]), Integer.parseInt(token[1]));
	}
	
	public void addScanNum(int scanNum)
	{
		if(specIndexList == null)
		{
			specIndexList = new ArrayList<Integer>();
			specIndexList.add(super.getFirst());
		}
		specIndexList.add(scanNum);
	}
	
	public ArrayList<Integer> getSpecIndexList()
	{
		return specIndexList;
	}
	
	public static ArrayList<SpecKey> getSpecKeyList(Iterator<Spectrum> itr, int minCharge, int maxCharge, ActivationMethod activationMethod)
	{
		if(activationMethod == ActivationMethod.FUSION)
			return getFusedSpecKeyList(itr, minCharge, maxCharge);
		
		ArrayList<SpecKey> specKeyList = new ArrayList<SpecKey>();
		
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			int specIndex = spec.getSpecIndex();
			int charge = spec.getCharge();
			
			if(activationMethod != null && spec.getActivationMethod() != null && spec.getActivationMethod() != activationMethod)
				continue;

			if(spec.size() < Constants.MIN_NUM_PEAKS_PER_SPECTRUM)
			{
//				System.out.println("Spectrum " + spec.getScanNum() + " has too few peaks (#Peaks: " + spec.size()+"): ignored.");
				continue;
			}
			if(charge == 0)
			{
				for(int c=minCharge; c<=maxCharge; c++)
					specKeyList.add(new SpecKey(specIndex, c));
			}
			else if(charge > 0)
			{
				specKeyList.add(new SpecKey(specIndex, charge));
			}
		}
		return specKeyList;
	}
	
	public static ArrayList<SpecKey> getFusedSpecKeyList(Iterator<Spectrum> itr, int minCharge, int maxCharge)
	{
		ArrayList<SpecKey> specKeyList = new ArrayList<SpecKey>();
		
		int prevSpecIndex = Integer.MIN_VALUE;
		int prevCharge = Integer.MIN_VALUE;
		float previousPrecursorMz = Float.MIN_VALUE;
		
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			int specIndex = spec.getSpecIndex();
			int charge = spec.getCharge();
			float precursorMz = spec.getPrecursorPeak().getMz();
			if(spec.getActivationMethod() == null)
			{
				System.out.println("Error: activation method is not available: Scan=" + spec.getSpecIndex()+", PrecursorMz=" + precursorMz);
				System.exit(-1);
			}
			
			if(specIndex == prevSpecIndex+1 && charge == prevCharge && precursorMz == previousPrecursorMz)
			{
				if(charge == 0)
				{
					for(int i=0; i<maxCharge-minCharge; i++)
						specKeyList.get(specKeyList.size()-1-i).addScanNum(specIndex);
				}
				else if(charge > 0)
				{
					specKeyList.get(specKeyList.size()-1).addScanNum(specIndex);
				}
			}
			else
			{
				if(charge == 0)
				{
					for(int c=minCharge; c<=maxCharge; c++)
						specKeyList.add(new SpecKey(specIndex, c));
				}
				else if(charge > 0)
				{
					specKeyList.add(new SpecKey(specIndex, charge));
				}
			}
			prevSpecIndex = specIndex;
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
			if(specKey.getSpecIndexList() == null || specKey.getSpecIndexList().size() != 2)
				System.out.println(specKey.getSpecKeyString()+"\t"+specKey.getSpecIndexList());
		}
		System.out.println("Size: " + list.size());
	}
	
}
