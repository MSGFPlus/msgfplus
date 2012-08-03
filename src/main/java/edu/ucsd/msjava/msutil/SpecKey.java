package edu.ucsd.msjava.msutil;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import edu.ucsd.msjava.parser.MzXMLSpectraIterator;
import edu.ucsd.msjava.sequences.Constants;


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
	
	public void addSpecIndex(int scanNum)
	{
		if(specIndexList == null)
		{
			specIndexList = new ArrayList<Integer>();
//			specIndexList.add(super.getFirst());
		}
		specIndexList.add(scanNum);
	}
	
	public ArrayList<Integer> getSpecIndexList()
	{
		return specIndexList;
	}
	
	public static ArrayList<SpecKey> getSpecKeyList(Iterator<Spectrum> itr, int startSpecIndex, int endSpecIndex, int minCharge, int maxCharge, ActivationMethod activationMethod)
	{
		if(activationMethod == ActivationMethod.FUSION)
			return getFusedSpecKeyList(itr, startSpecIndex, endSpecIndex, minCharge, maxCharge);
		
		ArrayList<SpecKey> specKeyList = new ArrayList<SpecKey>();
		
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			int specIndex = spec.getSpecIndex();
			
			if(specIndex < startSpecIndex)
				continue;
			if(specIndex >= endSpecIndex)
				continue;
			
			int charge = spec.getCharge();
			
			if(activationMethod != ActivationMethod.ASWRITTEN && spec.getActivationMethod() != null && spec.getActivationMethod() != activationMethod)
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
	
	public static ArrayList<SpecKey> getFusedSpecKeyList(Iterator<Spectrum> itr, int startSpecIndex, int endSpecIndex, int minCharge, int maxCharge)
	{
		HashMap<Peak, ArrayList<Integer>> precursorSpecIndexMap = new HashMap<Peak, ArrayList<Integer>>();
		
		while(itr.hasNext())
		{
			Spectrum spec = itr.next();
			int specIndex = spec.getSpecIndex();
			if(specIndex < startSpecIndex || specIndex >= endSpecIndex)
				continue;
			Peak precursor = spec.getPrecursorPeak();
			if(spec.getActivationMethod() == null)
			{
				System.out.println("Error: activation method is not available: Scan=" + spec.getSpecIndex()+", PrecursorMz=" + spec.getPrecursorPeak().getMz());
				System.exit(-1);
			}
			
			ArrayList<Integer> list = precursorSpecIndexMap.get(precursor);
			if(list == null)
			{
				list = new ArrayList<Integer>();
				precursorSpecIndexMap.put(precursor, list);
			}
			list.add(specIndex);
		}
		
		Iterator<Entry<Peak, ArrayList<Integer>>> mapItr = precursorSpecIndexMap.entrySet().iterator();
		ArrayList<SpecKey> specKeyList = new ArrayList<SpecKey>();
		while(mapItr.hasNext())
		{
			Entry<Peak, ArrayList<Integer>> entry = mapItr.next();
			Peak precursor = entry.getKey();
			ArrayList<Integer> list = entry.getValue();
			Collections.sort(list);
			
			int charge = precursor.getCharge();
			if(charge == 0)
			{
				for(int c=minCharge; c<=maxCharge; c++)
				{
					SpecKey specKey = new SpecKey(list.get(0), c);
					for(int specIndex : list)
						specKey.addSpecIndex(specIndex);
					specKeyList.add(specKey);
				}
			}
			else if(charge > 0)
			{
				SpecKey specKey = new SpecKey(list.get(0), charge);
				for(int specIndex : list)
					specKey.addSpecIndex(specIndex);
				specKeyList.add(specKey);
			}
			else
			{
				System.out.println("Error: negative precursor charge: " + precursor);
				System.exit(-1);
			}
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
		ArrayList<SpecKey> list = SpecKey.getFusedSpecKeyList(itr, 0, Integer.MAX_VALUE, minCharge, maxCharge);
		for(SpecKey specKey : list)
		{
			if(specKey.getSpecIndexList() == null || specKey.getSpecIndexList().size() != 2)
				System.out.println(specKey.getSpecKeyString()+"\t"+specKey.getSpecIndexList());
		}
		System.out.println("Size: " + list.size());
	}
	
//	public static ArrayList<SpecKey> getFusedSpecKeyListOld(Iterator<Spectrum> itr, int minCharge, int maxCharge)
//	{
//		ArrayList<SpecKey> specKeyList = new ArrayList<SpecKey>();
//		
//		int prevSpecIndex = Integer.MIN_VALUE;
//		int prevCharge = Integer.MIN_VALUE;
//		float previousPrecursorMz = Float.MIN_VALUE;
//		
//		while(itr.hasNext())
//		{
//			Spectrum spec = itr.next();
//			int specIndex = spec.getSpecIndex();
//			int charge = spec.getCharge();
//			float precursorMz = spec.getPrecursorPeak().getMz();
//			if(spec.getActivationMethod() == null)
//			{
//				System.out.println("Error: activation method is not available: Scan=" + spec.getSpecIndex()+", PrecursorMz=" + precursorMz);
//				System.exit(-1);
//			}
//			
//			if(specIndex == prevSpecIndex+1 && charge == prevCharge && precursorMz == previousPrecursorMz)
//			{
//				if(charge == 0)
//				{
//					for(int i=0; i<=maxCharge-minCharge; i++)
//						specKeyList.get(specKeyList.size()-1-i).addSpecIndex(specIndex);
//				}
//				else if(charge > 0)
//				{
//					specKeyList.get(specKeyList.size()-1).addSpecIndex(specIndex);
//				}
//			}
//			else
//			{
//				if(charge == 0)
//				{
//					for(int c=minCharge; c<=maxCharge; c++)
//					{
//						SpecKey specKey = new SpecKey(specIndex, c);
//						specKey.addSpecIndex(specIndex);
//						specKeyList.add(specKey);
//					}
//				}
//				else if(charge > 0)
//				{
//					SpecKey specKey = new SpecKey(specIndex, charge);
//					specKey.addSpecIndex(specIndex);
//					specKeyList.add(specKey);
//				}
//			}
//			prevSpecIndex = specIndex;
//			prevCharge = charge;
//			previousPrecursorMz = precursorMz;
//		}
//		return specKeyList;
//	}
	
}
