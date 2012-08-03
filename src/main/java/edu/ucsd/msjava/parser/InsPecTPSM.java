package edu.ucsd.msjava.parser;

import java.util.ArrayList;
import java.util.HashSet;

import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Annotation;
import edu.ucsd.msjava.msutil.Modification;
import edu.ucsd.msjava.msutil.ModifiedAminoAcid;
import edu.ucsd.msjava.msutil.Peptide;

public class InsPecTPSM extends PSM {
	private String insPecTString;
	private long specFilePos;
	private AminoAcid precedingAA;
	private AminoAcid succeedingAA;
	private ArrayList<Integer> scanNumList;
	
	public AminoAcid getPrecedingAA()	{ return precedingAA; }
	public AminoAcid getSucceedingAA()	{ return succeedingAA; }
	public Annotation getAnnotation()
	{
		return new Annotation(precedingAA, super.getPeptide(), succeedingAA);
	}
	
	public String getInsPecTString()	{ return insPecTString; }
	public long getSpecFilePos()	{ return specFilePos; }
	public ArrayList<Integer> getScanNumList() { return scanNumList; }
	
	public void setPrecedingAA(AminoAcid precedingAA)	{ this.precedingAA = precedingAA; }
	public void setSucceedingAA(AminoAcid succeedingAA)	{ this.succeedingAA = succeedingAA; }
	public void setSpecFilePos(long specFilePos) { this.specFilePos = specFilePos; }
	public void setInsPecTString(String insPecTString) { this.insPecTString = insPecTString; }
	public void setScanNumList(ArrayList<Integer> scanNumList) { this.scanNumList = scanNumList; }
	public AminoAcidSet getAASet(AminoAcidSet baseAASet)	
	{ 
		Peptide peptide = getPeptide();
		if(peptide != null && peptide.isModified())	
		{
			ArrayList<AminoAcid> modAAList = new ArrayList<AminoAcid>();
			for(AminoAcid aa : peptide)
			{
				if(aa.isModified())	// modified residue
				{
					modAAList.add(aa);
				}
			}
			return AminoAcidSet.getAminoAcidSetFromModAAList(baseAASet, modAAList);
		}
		else 
			return baseAASet;
	}

}
