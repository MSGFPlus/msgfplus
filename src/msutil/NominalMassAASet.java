package msutil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class NominalMassAASet {
	private Integer[] nominalMassArr;
	private HashMap<Integer,Integer> nominalMassNumTable;
	private HashMap<Integer,Float> nominalMassProbTable;
	private HashMap<Integer,ArrayList<AminoAcid>> nominalMassAATable;

	public NominalMassAASet(ArrayList<AminoAcid> aaList)
	{
		nominalMassNumTable = new HashMap<Integer, Integer>();
		nominalMassProbTable = new HashMap<Integer, Float>();
		nominalMassAATable = new HashMap<Integer, ArrayList<AminoAcid>>();
		for(AminoAcid aa : aaList)
		{
			float aaProb = aa.getProbability();
			int nominalMass = aa.getNominalMass();
			Integer nNum = nominalMassNumTable.get(nominalMass);
			if(nNum == null)
			{
				nominalMassNumTable.put(nominalMass, 1);
				nominalMassProbTable.put(nominalMass, aaProb);
				ArrayList<AminoAcid> newAAList = new ArrayList<AminoAcid>();
				newAAList.add(aa);
				nominalMassAATable.put(nominalMass, newAAList);
			}
			else
			{
				nominalMassNumTable.put(nominalMass, nNum+1);
				nominalMassProbTable.put(nominalMass, nominalMassProbTable.get(nominalMass)+aaProb);
				ArrayList<AminoAcid> existingAAList = nominalMassAATable.get(nominalMass);
				existingAAList.add(aa);
			}
		}		
		
		nominalMassArr = nominalMassNumTable.keySet().toArray(new Integer[0]);
		Arrays.sort(nominalMassArr);
	}
	
	/**
	 * Returns the array of distinct nominal masses
	 * @return the array of distinct nominal masses
	 */
	public Integer[] getDistinctNominalMasses() { return nominalMassArr; }
	
	/**
	 * Returns the number of amino acids with the given nominal mass
	 * @param nominalMass the nominal mass
	 * @return the number of amino acids with the given nominal mass
	 */
	public int getNumAAWithNominalMass(int nominalMass) 
	{ 
		Integer num = nominalMassNumTable.get(nominalMass);
		if(num == null)
			return 0;
		else
			return num;
	}

	/**
	 * Returns the sum of probabilities of amino acids with the given nominal mass
	 * @param nominalMass the nominal mass
	 * @return the sum of probabilities of amino acids with the given nominal mass
	 */
	public float getSumProbWithNominalMass(int nominalMass)
	{
		Float prob = nominalMassProbTable.get(nominalMass);
		if(prob == null)
			return 0;
		else
			return prob;
	}

	/**
	 * Returns the list of amino acids with the given nominal mass
	 * @param nominalMass the nominal mass
	 * @return the list of amino acids with the given nominal mass
	 */
	public ArrayList<AminoAcid> getAAListWithNominalMass(int nominalMass)
	{
		return nominalMassAATable.get(nominalMass);
	}
	
}
