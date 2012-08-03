package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.msutil.AminoAcidSet;

public class CandidatePeptideGridConsideringMetCleavage extends CandidatePeptideGrid {

	private final CandidatePeptideGrid candidatePepGridMetCleaved;		// For peptides with Met cleaved
	boolean isProteinNTermWithHeadingMet = false;
	
	public CandidatePeptideGridConsideringMetCleavage(AminoAcidSet aaSet, int maxPeptideLength) 
	{
		super(aaSet, maxPeptideLength);
		candidatePepGridMetCleaved = new CandidatePeptideGrid(aaSet, maxPeptideLength);
	}

	@Override
	public boolean addProtNTermResidue(char residue)
	{
		if(residue == 'M')
			isProteinNTermWithHeadingMet = true;
		else
			isProteinNTermWithHeadingMet = false;
		return super.addProtNTermResidue(residue);
	}

	@Override
	public boolean addNTermResidue(char residue)
	{
		isProteinNTermWithHeadingMet = false;
		return super.addNTermResidue(residue);
	}
	
	@Override
	public boolean addResidue(int length, char residue)
	{
		if(!super.addResidue(length, residue))
			return false;
		
		if(isProteinNTermWithHeadingMet)
		{
			if(length == 2)		// Second aa after M (e.g. _.M'G')
				return candidatePepGridMetCleaved.addProtNTermResidue(residue);
			else
				return candidatePepGridMetCleaved.addResidue(length-1, residue);
		}
		else
			return true;
	}
	
	@Override
	public boolean addProtCTermResidue(int length, char residue)
	{
		if(!super.addProtCTermResidue(length, residue))
			return false;
		
		if(isProteinNTermWithHeadingMet)
		{
			return candidatePepGridMetCleaved.addProtCTermResidue(length-1, residue);
		}
		else
			return true;
	}
	
	@Override
	public boolean addCTermResidue(int length, char residue)
	{
		if(!super.addCTermResidue(length, residue))
			return false;
		
		if(isProteinNTermWithHeadingMet)
		{
			return candidatePepGridMetCleaved.addCTermResidue(length-1, residue);
		}
		else
			return true;
	}	
	
	@Override
	public int size()
	{
		if(!isProteinNTermWithHeadingMet)
			return super.size();
		else
			return super.size()+candidatePepGridMetCleaved.size();
	}
	
	@Override
	public boolean isNTermMetCleaved(int index)
	{
		int sizeNormPep = super.size(); 
		if(index < sizeNormPep)
			return false;
		else
			return true;
	}
	
	@Override
	public int[] getNominalPRMGrid(int index)
	{
		if(!isProteinNTermWithHeadingMet)
			return super.getNominalPRMGrid(index);
		int sizeNormPep = super.size(); 
		if(index < sizeNormPep)
			return super.getNominalPRMGrid(index);
		else
			return candidatePepGridMetCleaved.getNominalPRMGrid(index-sizeNormPep);
	}
	
	@Override
	public double[] getPRMGrid(int index)
	{
		if(!isProteinNTermWithHeadingMet)
			return super.getPRMGrid(index);
		int sizeNormPep = super.size(); 
		if(index < sizeNormPep)
			return super.getPRMGrid(index);
		else
			return candidatePepGridMetCleaved.getPRMGrid(index-sizeNormPep);
	}
	
	@Override
	public float getPeptideMass(int index)
	{
		if(!isProteinNTermWithHeadingMet)
			return super.getPeptideMass(index); 
		int sizeNormPep = super.size(); 
		if(index < sizeNormPep)
			return super.getPeptideMass(index);
		else
			return candidatePepGridMetCleaved.getPeptideMass(index-sizeNormPep);
	}

	@Override
	public int getNominalPeptideMass(int index)
	{
		if(!isProteinNTermWithHeadingMet)
			return super.getNominalPeptideMass(index);
		int sizeNormPep = super.size(); 
		if(index < sizeNormPep)
			return super.getNominalPeptideMass(index);
		else
			return candidatePepGridMetCleaved.getNominalPeptideMass(index-sizeNormPep);
	}
	
	@Override
	public String getPeptideSeq(int index)
	{
		if(!isProteinNTermWithHeadingMet)
			return super.getPeptideSeq(index);
		int sizeNormPep = super.size(); 
		if(index < sizeNormPep)
			return super.getPeptideSeq(index);
		else
			return candidatePepGridMetCleaved.getPeptideSeq(index-sizeNormPep);
	}
	
	@Override
	public int getNumMods(int index)
	{
		if(!isProteinNTermWithHeadingMet)
			return super.getNumMods(index);
		int sizeNormPep = super.size(); 
		if(index < sizeNormPep)
			return super.getNumMods(index);
		else
			return candidatePepGridMetCleaved.getNumMods(index-sizeNormPep);
	}
}
