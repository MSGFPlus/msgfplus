package msutil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class CompositionAASet {
	private Composition[] compositionArr;
	private HashMap<Composition,Integer> compositionNumTable;
	private HashMap<Composition,Float> compositionProbTable;
	private HashMap<Composition,ArrayList<AminoAcid>> compositionAATable;

	public CompositionAASet(ArrayList<AminoAcid> aaList)
	{
		compositionNumTable = new HashMap<Composition, Integer>();
		compositionProbTable = new HashMap<Composition, Float>();
		compositionAATable = new HashMap<Composition, ArrayList<AminoAcid>>();

		for(AminoAcid aa : aaList)
		{
			float aaProb = aa.getProbability();
			
			Composition comp = aa.getComposition();
			Integer cNum = compositionNumTable.get(comp);
			if(cNum == null)
			{
				compositionNumTable.put(comp, 1);
				compositionProbTable.put(comp, aaProb);
				ArrayList<AminoAcid> newAAList = new ArrayList<AminoAcid>();
				newAAList.add(aa);
				compositionAATable.put(comp, newAAList);
			}
			else
			{
				compositionNumTable.put(comp, cNum+1);
				compositionProbTable.put(comp, compositionProbTable.get(aa.getComposition())+aaProb);
				ArrayList<AminoAcid> existingAAList = compositionAATable.get(comp);
				existingAAList.add(aa);
			}
		}
		
		compositionArr = compositionNumTable.keySet().toArray(new Composition[0]);
		Arrays.sort(compositionArr);		
	}
	
	/**
	 * Returns the array of distinct compositions
	 * @return the array of distinct compositions
	 */
	public Composition[] getDistinctCompositions() { return compositionArr; }
	
	/**
	 * Returns the number of amino acids with the given composition
	 * @param comp the composition
	 * @return the number of amino acids with the given composition
	 */
	public int getNumAAWithComposition(Composition comp) 
	{ 
		Integer num = compositionNumTable.get(comp);
		if(num == null)
			return 0;
		else
			return num;
	}

	/**
	 * Returns the sum of probabilities of amino acids with the given composition
	 * @param composition the composition
	 * @return the sum of probabilities of amino acids with the given composition
	 */
	public float getSumProbWithComposition(Composition comp)
	{
		Float prob = compositionProbTable.get(comp);
		if(prob == null)
			return 0;
		else
			return prob;
	}
	
	
	/**
	 * Returns the list of amino acids with the given composition
	 * @param composition the composition
	 * @return the list of amino acids with the given composition
	 */
	public ArrayList<AminoAcid> getIndexListWithComposition(Composition comp)
	{
		return compositionAATable.get(comp);
	}	
}
