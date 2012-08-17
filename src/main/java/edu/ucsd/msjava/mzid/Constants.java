package edu.ucsd.msjava.mzid;

import edu.ucsd.msjava.ui.MSGFPlus;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.Param;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;

public class Constants {
    static String analysisSoftID = "ID_software";
    static String providerID = "ID_provider";
    static String siListID = "SI_LIST_1";
    static String sirID = "SIR_";
    static String siiID = "SII_";
    static String spectraDataID = "SID_1";
    static String psiCvID = "PSI-MS";
    static String siProtocolID = "SearchProtocol_1";
    static String searchDBID = "SearchDB_1";
    static String pepEvidenceListID = "PepEvidList_1";
    static String specIdentID = "SpecIdent_1";
    static String unimodID = "UNIMOD";
    static String unitCvID = "UO";
    static String measureMzID = "Measure_MZ";
    static String measureIntID = "Measure_Int";
    static String measureErrorID = "Measure_Error";
    static String sourceFileID = "SourceFile_1";
	static final String UNIMOD_RESOURCE_PATH = "unimod.obo";

    static Cv psiCV;
    static Cv unimodCV;
    static Cv unitCV;

    static AnalysisSoftware msgfPlus;

//    static Person docOwner;
//    static Organization org;
//    static Affiliation aff;
    
    static {
		psiCV = new Cv();
		psiCV.setUri("http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");
		psiCV.setId(psiCvID);
		psiCV.setVersion("3.30.0");
		psiCV.setFullName("PSI-MS");

		unimodCV = new Cv();
		unimodCV.setUri("http://www.unimod.org/obo/unimod.obo");
		unimodCV.setId(unimodID);
		unimodCV.setFullName("UNIMOD");

		unitCV = new Cv();
		unitCV.setUri("http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo");
		unitCV.setId(unitCvID);
		unitCV.setFullName("UNIT-ONTOLOGY");
		
		msgfPlus = new AnalysisSoftware();
		msgfPlus.setName("MS-GF+");
		Param tempParam = new Param();
		tempParam.setParam(makeCvParam("MS:1002048","MS-GF+"));
//		tempParam.setParam(makeCvParam("MS:1001475","OMSSA"));
		msgfPlus.setSoftwareName(tempParam);
		msgfPlus.setId(analysisSoftID);
		msgfPlus.setVersion(MSGFPlus.VERSION);
		
//		docOwner = new Person();
//		docOwner.setId("PERSON_DOC_OWNER");
//		docOwner.setFirstName("Sangtae");
//		docOwner.setLastName("Kim");
//
//		org = new Organization();
//		org.setId("ORG_DOC_OWNER");
//		org.setName("UCSD");
//
//		Affiliation aff = new Affiliation();
//		aff.setOrganization(org);
//		docOwner.getAffiliation().add(aff);
		
    }
    
	/**
	 * Helper method to create and return a CvParam from accession, name and CV
	 *
	 * @return CvParam
	 */
	public static CvParam makeCvParam(String accession, String name){
		return makeCvParam(accession, name, psiCV);
	}
	
	/**
	 * Helper method to create and return a CvParam from accession, name and CV
	 *
	 * @return CvParam
	 */
	public static CvParam makeCvParam(String accession, String name, Cv cv){
		CvParam cvParam = new CvParam();
		cvParam.setAccession(accession);
		cvParam.setName(name);
		cvParam.setCv(cv);
		return cvParam;
	}

	/**
	 * Helper method to create and return a CvParam from accession, name, CV, unitAccession, unitName and unitCV 
	 *
	 * @return CvParam
	 */
	public static CvParam makeCvParam(String accession, String name, Cv cv, String unitAccession, String unitName, Cv alternateUnitCV){
		CvParam cvParam = makeCvParam(accession, name, cv);
		cvParam.setUnitAccession(unitAccession);
		cvParam.setUnitCv(alternateUnitCV);
		cvParam.setUnitName(unitName);
		return cvParam;
	}
	
	/**
	 * Helper method to create and return an UserParam
	 *
	 * @return UserParam
	 */
	public static UserParam makeUserParam(String name){
		UserParam userParam = new UserParam();
		userParam.setName(name);
		return userParam;
	}

	/**
	 * Helper method to create and return an UserParam
	 *
	 * @return UserParam
	 */
	public static UserParam makeUserParam(String name, String value){
		UserParam userParam = new UserParam();
		userParam.setName(name);
		userParam.setValue(value);
		return userParam;
	}
	
    /**
     * Helper method to setup a CvParam with CVRef, with either Daltons or ppm as units
     *
     */

    public static CvParam getCvParamWithMassUnits(Boolean isDaltonUnit){
        CvParam cvParam = new CvParam();

         //<cvParam accession="MS:1001413" name="search tolerance minus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
        cvParam.setCv(psiCV);
        cvParam.setUnitCv(unitCV);

        if(isDaltonUnit){
            cvParam.setUnitAccession("UO:0000221");
            cvParam.setUnitName("dalton");
        }
        else{
            cvParam.setUnitAccession("UO:0000169");
            cvParam.setUnitName("parts per million");
        }
        return cvParam;
    }
	
}
