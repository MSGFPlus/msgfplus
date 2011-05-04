package msdbsearch;

import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.ModifiedAminoAcid;

public class CandidatePeptideGrid {
	private static final int MAX_NUM_VARIATIONS_PER_PEPTIDE = 128;
	
	private final AminoAcidSet aaSet;
	private final int maxPeptideLength;
	private final int numMaxMods;
	
	private int[][] nominalPRM;
	private double[][] prm;
	private int[][] numMods;
	private StringBuffer[] peptide;
	
	// anywhere aa (including modified aa)
	private int[][] aaNominalMass; // residue -> mass list
	private double[][] aaMass;
	private char[][] aaResidue;
	
	// aa with variable modifications
	private int[][] nTermAANominalMass; // residue -> mass list
	private double[][] nTermAAMass;
	private char[][] nTermAAResidue;
	
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
		//TODO: deal with non aa-specific N-term mods 
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
	
	public boolean addResidue(char residue)
	{
		return addResidue(length+1, residue);
	}
	
	public int[][] getNominalPRMGrid()
	{
		return this.nominalPRM;
	}
	
	public double[][] getPRMGrid()
	{
		return this.prm;
	}
	
	public int size()
	{
		return size[length];
	}
	
	public float getPeptideMass(int index)
	{
		return (float)prm[index][length];
	}
	
	public String getPeptideSeq(int index)
	{
		return peptide[index].toString();
	}
	
	// if residue is not a standard residue, return false
	public boolean addResidue(int length, char residue)
	{
		double[] aaMassArr;
		if(length != 1)
			aaMassArr = aaMass[residue];
		else	// N-term
			aaMassArr = nTermAAMass[residue];
		if(aaMassArr == null || length > maxPeptideLength)
			return false;
		
		int[] aaNominalMassArr;
		char[] aaResidueArr;
		if(length != 1)
		{
			aaNominalMassArr = aaNominalMass[residue];
			aaResidueArr = aaResidue[residue];
		}
		else	// N-term
		{
			aaNominalMassArr = nTermAANominalMass[residue];
			aaResidueArr = nTermAAResidue[residue];
		}
		
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
		aaNominalMass = new int[128][];
		aaMass = new double[128][];
		aaResidue = new char[128][];
		for(AminoAcid aa : AminoAcidSet.getStandardAminoAcidSet())
		{
			char residue = aa.getResidue();
			AminoAcid[] aaArr = aaSet.getAminoAcids(residue);
			aaNominalMass[residue] = new int[aaArr.length];
			aaMass[residue] = new double[aaArr.length];
			aaResidue[residue] = new char[aaArr.length];
			for(int i=0; i<aaArr.length; i++)
			{
				aaNominalMass[residue][i] = aaArr[i].getNominalMass();
				aaMass[residue][i] = aaArr[i].getAccurateMass();
				aaResidue[residue][i] = aaArr[i].getResidue();
			}
		}
		
		// N-term
		nTermAANominalMass = new int[128][];
		nTermAAMass = new double[128][];
		nTermAAResidue = new char[128][];
		for(AminoAcid aa : AminoAcidSet.getStandardAminoAcidSet())
		{
			char residue = aa.getResidue();
			AminoAcid[] aaArr = aaSet.getNTermAminoAcids(residue);
			nTermAANominalMass[residue] = new int[aaArr.length];
			nTermAAMass[residue] = new double[aaArr.length];
			nTermAAResidue[residue] = new char[aaArr.length];
			for(int i=0; i<aaArr.length; i++)
			{
				nTermAANominalMass[residue][i] = aaArr[i].getNominalMass();
				nTermAAMass[residue][i] = aaArr[i].getAccurateMass();
				nTermAAResidue[residue][i] = aaArr[i].getResidue();
			}
		}
	}	
}
