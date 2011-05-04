package parser;

import java.util.ArrayList;

import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Peptide;


public class MSGappedDictionaryPSM  extends PSM{
	private AminoAcid precedingAA;
	private AminoAcid succeedingAA;
	private AminoAcidSet aaSet;
	private float parentMassError;
	//private boolean isPeptideModified = false; // remove next. PSM will take care of this
	
	public MSGappedDictionaryPSM aaSet(AminoAcidSet aaSet)	{ this.aaSet = aaSet; return this;}
	
	public MSGappedDictionaryPSM peptide(String peptideStr) {
		ArrayList<AminoAcid> aaList = new ArrayList<AminoAcid>();
		for(int i=0; i<peptideStr.length();i++)
			aaList.add(aaSet.getAminoAcid(peptideStr.charAt(i)));
		this.peptide(new Peptide(aaList));
		return this;
	}
	
	public boolean isPeptideModified(){
		for(int i=0; i<this.getPeptideStr().length(); i++)
			if(Character.isLowerCase(this.getPeptideStr().charAt(i))) return true;
		return false;
	}
	
	public AminoAcid getPrecedingAA()	{ return precedingAA; }
	public AminoAcid getSucceedingAA()	{ return succeedingAA; }
	public AminoAcidSet getAASet()	{ return aaSet; }
	public float getParentMassError() {return parentMassError;}
	
	public void setPrecedingAA(AminoAcid precedingAA)	{ this.precedingAA = precedingAA; }
	public void setSucceedingAA(AminoAcid succeedingAA)	{ this.succeedingAA = succeedingAA; }
	public void setParentMassError(float error) {this.parentMassError = error;}
}
