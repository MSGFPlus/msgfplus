package edu.ucsd.msjava.ipa;

import edu.ucsd.msjava.mzid.UnimodComposition;

public class PSM {
	PSM(String resultStr)
	{
		String[] token = resultStr.split("\t");
		this.scanNum = Integer.parseInt(token[2]);
		this.precursorMz = Float.parseFloat(token[3]);
		this.charge = Integer.parseInt(token[7]);
		this.peptide = token[8];
		this.composition = new UnimodComposition();
		composition.add(token[9]);
		this.specEValue = Float.parseFloat(token[13]);
	}

	int getScanNum() {
		return scanNum;
	}
	float getPrecursorMz() {
		return precursorMz;
	}
	int getCharge() {
		return charge;
	}
	String getPeptide() {
		return peptide;
	}
	UnimodComposition getComposition() {
		return composition;
	}
	float getSpecEValue() {
		return specEValue;
	}
	
	private int scanNum;
	private float precursorMz;
	private int charge;
	private String peptide;
	private UnimodComposition composition;
	private float specEValue;
}
