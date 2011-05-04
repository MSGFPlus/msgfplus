package msgf;

import java.util.HashSet;
import java.util.Hashtable;

public class ReachableNode {
	public static Hashtable<Integer, HashSet<Integer>> possibleNominalMassTable = new Hashtable<Integer, HashSet<Integer>>();
	
	public static HashSet<Integer> getPossibleNominalMassSet(int nominalPeptideMass)
	{
		HashSet<Integer> set = possibleNominalMassTable.get(nominalPeptideMass);
		return set;
	}
}
