package edu.ucsd.msjava.msutil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import edu.ucsd.msjava.msgf.IntMassFactory;

public class IntMassAASet {
	private Integer[] intMassArr;
	private HashMap<Integer,Integer> intMassNumTable;
	private HashMap<Integer,Float> intMassProbTable;
	private HashMap<Integer,ArrayList<AminoAcid>> intMassAATable;

	public IntMassAASet(IntMassFactory factory, ArrayList<AminoAcid> aaList)
	{
		intMassNumTable = new HashMap<Integer, Integer>();
		intMassProbTable = new HashMap<Integer, Float>();
		intMassAATable = new HashMap<Integer, ArrayList<AminoAcid>>();
		for(AminoAcid aa : aaList)
		{
			float aaProb = aa.getProbability();
			int intMass = factory.getMassIndex(aa.getMass());
			Integer nNum = intMassNumTable.get(intMass);
			if(nNum == null)
			{
				intMassNumTable.put(intMass, 1);
				intMassProbTable.put(intMass, aaProb);
				ArrayList<AminoAcid> newAAList = new ArrayList<AminoAcid>();
				newAAList.add(aa);
				intMassAATable.put(intMass, newAAList);
			}
			else
			{
				intMassNumTable.put(intMass, nNum+1);
				intMassProbTable.put(intMass, intMassProbTable.get(intMass)+aaProb);
				ArrayList<AminoAcid> existingAAList = intMassAATable.get(intMass);
				existingAAList.add(aa);
			}
		}		
		
		intMassArr = intMassNumTable.keySet().toArray(new Integer[0]);
		Arrays.sort(intMassArr);
	}
	
	/**
	 * Returns the array of distinct integer masses
	 * @return the array of distinct integer masses
	 */
	public Integer[] getDistinctIntMasses() { return intMassArr; }
	
	/**
	 * Returns the number of amino acids with the given integer mass
	 * @param intMass the integer mass
	 * @return the number of amino acids with the given integer mass
	 */
	public int getNumAAWithIntMass(int intMass) 
	{ 
		Integer num = intMassNumTable.get(intMass);
		if(num == null)
			return 0;
		else
			return num;
	}

	/**
	 * Returns the sum of probabilities of amino acids with the given integer mass
	 * @param intMass the integer mass
	 * @return the sum of probabilities of amino acids with the given integer mass
	 */
	public float getSumProbWithIntMass(int intMass)
	{
		Float prob = intMassProbTable.get(intMass);
		if(prob == null)
			return 0;
		else
			return prob;
	}

	/**
	 * Returns the list of amino acids with the given integer mass
	 * @param intMass the integer mass
	 * @return the list of amino acids with the given integer mass
	 */
	public ArrayList<AminoAcid> getAAListWithIntMass(int intMass)
	{
		return intMassAATable.get(intMass);
	}
	
}
