package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Modification.Location;

public class CandidatePeptideGrid {
	private static final int MAX_NUM_VARIATIONS_PER_PEPTIDE = 128;
	private static final int STANDARD_RESIDUE_MAX_RESIDUE = 128;
	
	private final AminoAcidSet aaSet;
	private final int maxPeptideLength;
	private final int numMaxMods;
	
	private int[][] nominalPRM;
	private double[][] prm;
	private int[][] numMods;
	private StringBuffer[] peptide;
	
	// caching amino acid set for fast search
	
	// anywhere aa (including modified aa)
	private int[][] aaNominalMass; // residue -> mass list
	private double[][] aaMass;
	private char[][] aaResidue;
	
	// N-term aa set
	private int[][] nTermAANominalMass; // residue -> mass list
	private double[][] nTermAAMass;
	private char[][] nTermAAResidue;

	// C-term aa set
	private int[][] cTermAANominalMass; // residue -> mass list
	private double[][] cTermAAMass;
	private char[][] cTermAAResidue;

	// Protein N-term aa set
	private int[][] protNTermAANominalMass; // residue -> mass list
	private double[][] protNTermAAMass;
	private char[][] protNTermAAResidue;

	// Protein C-term aa set
	private int[][] protCTermAANominalMass; // residue -> mass list
	private double[][] protCTermAAMass;
	private char[][] protCTermAAResidue;
	
	// Protein N-term Met cleavage
	private int length;
	private int[] size;
	
	public CandidatePeptideGrid(AminoAcidSet aaSet, int maxPeptideLength)
	{
		this.numMaxMods = aaSet.getMaxNumberOfVariableModificationsPerPeptide();
		this.maxPeptideLength = maxPeptideLength;
		this.aaSet = aaSet;
		
		cacheAASet();
		
		nominalPRM = new int[MAX_NUM_VARIATIONS_PER_PEPTIDE][maxPeptideLength+1];
		prm = new double[MAX_NUM_VARIATIONS_PER_PEPTIDE][maxPeptideLength+1];
		numMods = new int[MAX_NUM_VARIATIONS_PER_PEPTIDE][maxPeptideLength+1];
		peptide = new StringBuffer[MAX_NUM_VARIATIONS_PER_PEPTIDE];
		size = new int[MAX_NUM_VARIATIONS_PER_PEPTIDE];
		
		initializeNTerm();
	}	
	
	private void initializeNTerm()
	{
		for(int i=0; i<MAX_NUM_VARIATIONS_PER_PEPTIDE; i++)
		{
			nominalPRM[i][0] = 0;
			prm[i][0] = 0.;
			numMods[i][0] = 0;
			peptide[i] = new StringBuffer();
		}
		size[0] = 1;
		length = 0;
	}
	
	public int[] getNominalPRMGrid(int index)
	{
		return this.nominalPRM[index];
	}
	
	public double[] getPRMGrid(int index)
	{
		return this.prm[index];
	}
	
	public int size()
	{
		return size[length];
	}
	
	public float getPeptideMass(int index)
	{
		return (float)prm[index][length];
	}

	public int getNominalPeptideMass(int index)
	{
		return nominalPRM[index][length];
	}
	
	public String getPeptideSeq(int index)
	{
		return peptide[index].toString();
	}
	
	public int getNumMods(int index)
	{
		return numMods[index][length];
	}
	
//	public boolean addResidue(char residue)
//	{
//		return addResidue(length+1, residue);
//	}
	
	// if residue is not a standard residue, return false
	public boolean addResidue(int length, char residue)
	{
		double[] aaMassArr = aaMass[residue];
		if(aaMassArr == null || length > maxPeptideLength)
			return false;
		
		int[] aaNominalMassArr = aaNominalMass[residue];
		char[] aaResidueArr = aaResidue[residue];
		
		return addResidue(aaMassArr, aaNominalMassArr, aaResidueArr, length);
	}
	
	public boolean addProtNTermResidue(char residue)
	{
		double[] aaMassArr = protNTermAAMass[residue];
		if(aaMassArr == null)
			return false;
		
		int[] aaNominalMassArr = protNTermAANominalMass[residue];
		char[] aaResidueArr = protNTermAAResidue[residue];
		
		return addResidue(aaMassArr, aaNominalMassArr, aaResidueArr, 1);
	}

	public boolean addNTermResidue(char residue)
	{
		double[] aaMassArr = nTermAAMass[residue];
		if(aaMassArr == null)
			return false;
		
		int[] aaNominalMassArr = nTermAANominalMass[residue];
		char[] aaResidueArr = nTermAAResidue[residue];
		
		return addResidue(aaMassArr, aaNominalMassArr, aaResidueArr, 1);
	}
	
	public boolean addProtCTermResidue(int length, char residue)
	{
		double[] aaMassArr = protCTermAAMass[residue];
		if(aaMassArr == null)
			return false;
		
		int[] aaNominalMassArr = protCTermAANominalMass[residue];
		char[] aaResidueArr = protCTermAAResidue[residue];
		
		return addResidue(aaMassArr, aaNominalMassArr, aaResidueArr, length);
	}
	
	public boolean addCTermResidue(int length, char residue)
	{
		double[] aaMassArr = cTermAAMass[residue];
		if(aaMassArr == null)
			return false;
		
		int[] aaNominalMassArr = cTermAANominalMass[residue];
		char[] aaResidueArr = cTermAAResidue[residue];
		
		return addResidue(aaMassArr, aaNominalMassArr, aaResidueArr, length);
	}
	
	public boolean isNTermMetCleaved(int index)	{ return false; }
	
	private boolean addResidue(double[] aaMassArr, int[] aaNominalMassArr, char[] aaResidueArr, int length)
	{
		int parentSize = size[length-1];
		for(int parentIndex=0; parentIndex<parentSize; parentIndex++)
		{
			nominalPRM[parentIndex][length] = nominalPRM[parentIndex][length-1] + aaNominalMassArr[0];
			prm[parentIndex][length] = prm[parentIndex][length-1] + aaMassArr[0];
			numMods[parentIndex][length] = numMods[parentIndex][length-1];
			peptide[parentIndex].setLength(length-1);
			peptide[parentIndex].append(aaResidueArr[0]);
		}
		size[length] = parentSize;
		
		// modified residue: copy prms up to length-1 into new array		
		if(aaMassArr.length > 1 && parentSize < MAX_NUM_VARIATIONS_PER_PEPTIDE)	
		{
			int newIndex = parentSize;
			for(int parentIndex=0; parentIndex<parentSize; parentIndex++)
			{
				int numModParent = numMods[parentIndex][length-1]; 
				if(numModParent < numMaxMods)
				{
					for(int j=1; j<aaMassArr.length; j++)
					{
						for(int k=1; k<length; k++)
						{
							nominalPRM[newIndex][k] = nominalPRM[parentIndex][k];
							prm[newIndex][k] = prm[parentIndex][k];
						}
						peptide[newIndex] = new StringBuffer(peptide[parentIndex].substring(0, length-1));
						nominalPRM[newIndex][length] = nominalPRM[newIndex][length-1] + aaNominalMassArr[j];
						prm[newIndex][length] = prm[newIndex][length-1] + aaMassArr[j];
						numMods[newIndex][length] = numModParent+1;
						peptide[newIndex].append(aaResidueArr[j]);
						newIndex++;
						if(newIndex >= MAX_NUM_VARIATIONS_PER_PEPTIDE)
							break;
					}
				}
				if(newIndex >= MAX_NUM_VARIATIONS_PER_PEPTIDE)
					break;
			}
			size[length] = newIndex;
		}
		
		this.length = length;
		return true;		
	}
	
	private void cacheAASet()
	{
		for(Location location : Location.values())
			cacheAASet(location);
	}
	
	private void cacheAASet(Location location)
	{
		int[][] stdResidue2NominalMasses = null;
		double[][] stdResidue2Masses = null;
		char[][] stdResidue2Residues = null;
		
		if(location == Location.Anywhere)
		{
			stdResidue2NominalMasses = aaNominalMass = new int[STANDARD_RESIDUE_MAX_RESIDUE][];
			stdResidue2Masses = aaMass = new double[STANDARD_RESIDUE_MAX_RESIDUE][];
			stdResidue2Residues = aaResidue = new char[STANDARD_RESIDUE_MAX_RESIDUE][];
		}
		else if(location == Location.N_Term)
		{
			stdResidue2NominalMasses = nTermAANominalMass = new int[STANDARD_RESIDUE_MAX_RESIDUE][];
			stdResidue2Masses = nTermAAMass = new double[STANDARD_RESIDUE_MAX_RESIDUE][];
			stdResidue2Residues = nTermAAResidue = new char[STANDARD_RESIDUE_MAX_RESIDUE][];
		}
		else if(location == Location.C_Term)
		{
			stdResidue2NominalMasses = cTermAANominalMass = new int[STANDARD_RESIDUE_MAX_RESIDUE][];
			stdResidue2Masses = cTermAAMass = new double[STANDARD_RESIDUE_MAX_RESIDUE][];
			stdResidue2Residues = cTermAAResidue = new char[STANDARD_RESIDUE_MAX_RESIDUE][];
		}
		else if(location == Location.Protein_N_Term)
		{
			stdResidue2NominalMasses = protNTermAANominalMass = new int[STANDARD_RESIDUE_MAX_RESIDUE][];
			stdResidue2Masses = protNTermAAMass = new double[STANDARD_RESIDUE_MAX_RESIDUE][];
			stdResidue2Residues = protNTermAAResidue = new char[STANDARD_RESIDUE_MAX_RESIDUE][];
		}
		else if(location == Location.Protein_C_Term)
		{
			stdResidue2NominalMasses = protCTermAANominalMass = new int[STANDARD_RESIDUE_MAX_RESIDUE][];
			stdResidue2Masses = protCTermAAMass = new double[STANDARD_RESIDUE_MAX_RESIDUE][];
			stdResidue2Residues = protCTermAAResidue = new char[STANDARD_RESIDUE_MAX_RESIDUE][];
		}
		
		for(AminoAcid aa : AminoAcidSet.getStandardAminoAcidSet())
		{
			char residue = aa.getResidue();
			AminoAcid[] aaArr = aaSet.getAminoAcids(location, residue);
			stdResidue2NominalMasses[residue] = new int[aaArr.length];
			stdResidue2Masses[residue] = new double[aaArr.length];
			stdResidue2Residues[residue] = new char[aaArr.length];
			for(int i=0; i<aaArr.length; i++)
			{
				stdResidue2NominalMasses[residue][i] = aaArr[i].getNominalMass();
				stdResidue2Masses[residue][i] = aaArr[i].getAccurateMass();
				stdResidue2Residues[residue][i] = aaArr[i].getResidue();
			}
		}
	}
}
