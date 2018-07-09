package edu.ucsd.msjava.msutil;

/**
 * Stored information for a cvParam, for "copying" cvParams from mzML( and potentially other formats) to mzIdentML
 * @author Bryson Gibbons
 */
public class CvParamInfo {
    private String accession;
    private String name;
    private String value;
    private String unitAccession;
    private String unitName;
    private Boolean hasUnit;

    public CvParamInfo(String accession, String name, String value) {
        this.accession = accession;
        this.name = name;
        this.value = value;
    }

    public CvParamInfo(String accession, String name, String value, String unitAccession, String unitName) {
        this.accession = accession;
        this.name = name;
        this.value = value;
        this.hasUnit = true;
        this.unitAccession = unitAccession;
        this.unitName = unitName;
    }

    public String getAccession() {
        return this.accession;
    }

    public String getName() {
        return this.name;
    }

    public String getValue() {
        return this.value;
    }

    public Boolean getHasUnit() {
        return this.hasUnit;
    }

    public String getUnitAccession() {
        return this.unitAccession;
    }

    public String getUnitName() {
        return this.unitName;
    }
}
