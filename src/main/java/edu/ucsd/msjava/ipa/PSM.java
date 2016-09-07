package edu.ucsd.msjava.ipa;

import edu.ucsd.msjava.mzid.UnimodComposition;

public class PSM {
    PSM(String resultStr) {
        this.resultStr = resultStr;
        String[] token = resultStr.split("\t");
        this.scanNum = Integer.parseInt(token[2]);
        this.precursorMz = Float.parseFloat(token[5]);
        this.charge = Integer.parseInt(token[8]);
        this.peptide = token[9];
        this.composition = new UnimodComposition();
        composition.add(token[10]);
        this.specEValue = Float.parseFloat(token[14]);
        this.eValue = Float.parseFloat(token[15]);
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

    float getEValue() {
        return eValue;
    }

    String getResultString() {
        return resultStr;
    }

    private String resultStr;
    private int scanNum;
    private float precursorMz;
    private int charge;
    private String peptide;
    private UnimodComposition composition;
    private float specEValue;
    private float eValue;
}
