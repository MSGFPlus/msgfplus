package parser;

import java.util.ArrayList;

import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Annotation;

public class InsPecTPSM extends PSM {
	private String insPecTString;
	private long specFilePos;
	private AminoAcid precedingAA;
	private AminoAcid succeedingAA;
	private AminoAcidSet aaSet;
	private ArrayList<Integer> scanNumList;
	
	public AminoAcid getPrecedingAA()	{ return precedingAA; }
	public AminoAcid getSucceedingAA()	{ return succeedingAA; }
	public Annotation getAnnotation()
	{
		return new Annotation(precedingAA, super.getPeptide(), succeedingAA);
	}
	
	public String getInsPecTString()	{ return insPecTString; }
	public long getSpecFilePos()	{ return specFilePos; }
	public AminoAcidSet getAASet()	{ return aaSet; }
	public ArrayList<Integer> getScanNumList() { return scanNumList; }
	
	public void setPrecedingAA(AminoAcid precedingAA)	{ this.precedingAA = precedingAA; }
	public void setSucceedingAA(AminoAcid succeedingAA)	{ this.succeedingAA = succeedingAA; }
	public void setAASet(AminoAcidSet aaSet)	{ this.aaSet = aaSet; }
	public void setSpecFilePos(long specFilePos) { this.specFilePos = specFilePos; }
	public void setInsPecTString(String insPecTString) { this.insPecTString = insPecTString; }
	public void setScanNumList(ArrayList<Integer> scanNumList) { this.scanNumList = scanNumList; }
}
