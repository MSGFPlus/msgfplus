package msdbsearch;

import java.util.ArrayList;

import msutil.AminoAcid;
import msutil.AminoAcidSet;

public class CandidatePeptideTree {
	private static final int MAX_NUM_VARIATIONS_PER_PEPTIDE = 128;
	
	private final int maxPeptideLength;
	private final AminoAcidSet aaSet;
	private final int numMaxMods;
	
	private int length;
	private Node[][] treeNodes;
	
	private int[][] nominalPRM;
	private double[][] prm;
	private StringBuffer[] peptide;
	
	public CandidatePeptideTree(AminoAcidSet aaSet, int maxPeptideLength)
	{
		this.aaSet = aaSet;
		this.numMaxMods = aaSet.getMaxNumberOfVariableModificationsPerPeptide();
		this.maxPeptideLength = maxPeptideLength;
		
		treeNodes = new Node[maxPeptideLength+1][];
		treeNodes[0] = new Node[1];
		treeNodes[0][0] = Node.ROOT;
		
		nominalPRM = new int[maxPeptideLength+1][MAX_NUM_VARIATIONS_PER_PEPTIDE];
		prm = new double[maxPeptideLength+1][MAX_NUM_VARIATIONS_PER_PEPTIDE];
		peptide = new StringBuffer[MAX_NUM_VARIATIONS_PER_PEPTIDE];
		
		length = 0;
	}	
	
	public Node[] getNodesOfLength(int length)
	{
		return treeNodes[length];
	}
	
	public boolean addResidue(char residue)
	{
		return addResidue(length+1, residue);
	}
	
	// if residue is not a standard residue, return false
	public boolean addResidue(int length, char residue)
	{
		AminoAcid[] aaArr = aaSet.getAminoAcids(residue);
		if(aaArr == null)
			return false;
		Node[] parents = treeNodes[length-1];
		ArrayList<Node> newNodes = new ArrayList<Node>();
		for(int i=0; i<parents.length; i++)
		{
			Node parent = parents[i];
			int numModeParent = parent.numMods;
			boolean addModifiedAA;
			if(numModeParent < numMaxMods)
				addModifiedAA = true;
			else
				addModifiedAA = false;
			for(int j=0; j<aaArr.length; j++)
			{
				AminoAcid aa = aaArr[j];
				if(addModifiedAA || !aa.isModified())
					newNodes.add(new Node(parent, aaArr[j]));
			}
			treeNodes[length] = newNodes.toArray(new Node[0]);
		}
		this.length = length;
		return true;
	}
	
	public ArrayList<CandidatePeptide> getCandidatePeptides()
	{
		ArrayList<CandidatePeptide> candidates = new ArrayList<CandidatePeptide>();
		for(Node terminalNode : treeNodes[length])
		{
			CandidatePeptide candidate = new CandidatePeptide(length);
			Node curNode = terminalNode;
			while(curNode != Node.ROOT)
			{
				candidate.add(curNode.aa);
				curNode = curNode.parent;
			}
		}
		return candidates;
	}
	
	public static class CandidatePeptide {
		private double[] srm;
		private int[] nominalSRM;
		int length = 0;
		public CandidatePeptide(int length)
		{
			srm = new double[length+1];
			nominalSRM = new int[length+1];
			srm[0] = 0.;
			nominalSRM[0] = 0;
			length = 0;
		}
		public void add(AminoAcid aa)
		{
			srm[length+1] = srm[length]+aa.getAccurateMass();
			nominalSRM[length+1] = nominalSRM[length]+aa.getNominalMass();
			length++;
		}
		
		public double[] getSRMArray()
		{
			return srm;
		}
		public int[] getNominalSRMArray()
		{
			return nominalSRM;
		}
		public float getPeptideMass()
		{
			return 0;
		}
	}
	
	private static class Node {
		Node parent;
		AminoAcid aa;
//		double prm;
//		int nominalPRM;
		int numMods;
		
		public static final Node ROOT = new Node();
		private Node() 
		{
			parent = null;
			aa = null;
//			prm = 0.;
//			nominalPRM = 0;
			numMods = 0;
		}
		public Node(Node parent, AminoAcid aa)
		{
			this.parent = parent;
			this.aa = aa;
//			this.prm = parent.prm + aa.getAccurateMass();
//			this.nominalPRM = parent.nominalPRM + aa.getNominalMass();
			if(aa.isModified())
				numMods = parent.numMods+1;
			else
				numMods = parent.numMods;
		}
	}
	
}
