package parser;

import java.io.PrintStream;
import java.util.List;

import uk.ac.ebi.jmzidml.model.mzidml.AbstractContact;
import uk.ac.ebi.jmzidml.model.mzidml.Affiliation;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.ContactRole;
import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DataCollection;
import uk.ac.ebi.jmzidml.model.mzidml.Organization;
import uk.ac.ebi.jmzidml.model.mzidml.Param;
import uk.ac.ebi.jmzidml.model.mzidml.Person;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;
import uk.ac.ebi.jmzidml.model.mzidml.Role;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;

public class MZIdentMLGen {
	private static String psiCvID = "PSI-MS";
	private static String unimodID = "UNIMOD";
	private static String unitCvID = "UO";
	private static String analysisSoftID = "ID_software";
   
    private MzIdentMLMarshaller m;
    private Cv psiCV;
    private Person docOwner;
    private Organization org;
    
    
    public MZIdentMLGen()
    {
    	m = new MzIdentMLMarshaller();
    	
        Cv psiCV = new Cv();
        psiCV.setUri("http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");
        psiCV.setId(psiCvID);
        psiCV.setVersion("2.25.0");
        psiCV.setFullName("PSI-MS");
        
    }

    public void writeResults(PrintStream out)
    {
        out.println(m.createXmlHeader());
        out.println(m.createMzIdentMLStartTag("12345"));
    	m.marshal(generateCvList(), out);
    	out.println();
    	m.marshal(generateAnalysisSoftwareList(), out);
    	out.println();
    	
    	Provider provider = generateProvider();
    	if(provider != null)
    	{
    		m.marshal(generateProvider(), out);
        	out.println();
    	}
    	AuditCollection auditCollection = generateAuditCollection();
    	if(auditCollection != null)
    	{
    		m.marshal(generateAuditCollection(), out);
        	out.println();
    	}
    	
    	m.marshal(generateSequenceCollection(), out);
    	out.println();

    	m.marshal(generateAnalysisCollection(), out);
    	out.println();

    	m.marshal(generateAnalysisProtocolCollection(), out);
    	out.println();

    	m.marshal(generateDataCollection(), out);
    	out.println();
    }
    
    public MZIdentMLGen setPersonalInfo(
    		String firstName,
    		String lastName,
    		String affiliationName
    		)
    {
        docOwner = new Person();
        docOwner.setId("PERSON_DOC_OWNER");
        docOwner.setFirstName(firstName);
        docOwner.setLastName(lastName);
        
        org = new Organization();
        org.setId("ORG_DOC_OWNER");
        org.setName(affiliationName);

        Affiliation aff = new Affiliation();
        aff.setOrganization(org);
        docOwner.getAffiliation().add(aff);
        
    	return this;
    }
    
    private CvList generateCvList()
    {
        CvList cvList = new CvList();
        List<Cv> localCvList = cvList.getCv();

        Cv unimodCV = new Cv();
        unimodCV.setUri("http://www.unimod.org/obo/unimod.obo");
        unimodCV.setId(unimodID);
        unimodCV.setFullName("UNIMOD");

        Cv unitCV = new Cv();
        unitCV.setUri("http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo");
        unitCV.setId(unitCvID);
        unitCV.setFullName("UNIT-ONTOLOGY");

        localCvList.add(psiCV);
        localCvList.add(unimodCV);
        localCvList.add(unitCV);
        
        return cvList;
    }
    
    private AnalysisSoftwareList generateAnalysisSoftwareList()
    {
        AnalysisSoftwareList analysisSoftwareList = new AnalysisSoftwareList();
        List<AnalysisSoftware> analysisSoftwares = analysisSoftwareList.getAnalysisSoftware();
        AnalysisSoftware analysisSoftware = new AnalysisSoftware();
        analysisSoftware.setName("MS-GF+");
        Param tempParam = new Param();
        tempParam.setParam(makeCvParam("MS:XXXXXXX","MS-GF+",psiCV));
        analysisSoftware.setSoftwareName(tempParam);
        analysisSoftware.setId(analysisSoftID);
        analysisSoftware.setVersion("1.0");
        analysisSoftwares.add(analysisSoftware);
        
        return analysisSoftwareList;    	
    }

    private Provider generateProvider()
    {
    	if(docOwner == null)
    		return null;
    	
    	Provider provider = new Provider();
        provider.setId("PROVIDER");

        ContactRole contactRole = new ContactRole();
        contactRole.setContact(docOwner);

        Role role = new Role();
        role.setCvParam(makeCvParam("MS:1001271","researcher",psiCV));
        contactRole.setRole(role);
        provider.setContactRole(contactRole);
        
        return provider;
    }
    
    
    private AuditCollection generateAuditCollection()
    {
    	if(docOwner == null || org == null)
    		return null;
    	
    	AuditCollection auditCollection = new AuditCollection();
        List<AbstractContact> contactList = auditCollection.getPersonOrOrganization();
        
       	contactList.add(docOwner);
       	contactList.add(org);

        return auditCollection;
    }
    
    private SequenceCollection generateSequenceCollection()
    {
    	SequenceCollection sequenceCollection = new SequenceCollection();
    	return sequenceCollection;
    }

    private AnalysisCollection generateAnalysisCollection()
    {
    	AnalysisCollection analysisCollection = new AnalysisCollection();
    	return analysisCollection;
    }

    private AnalysisProtocolCollection generateAnalysisProtocolCollection()
    {
    	AnalysisProtocolCollection analysisProtocolCollection = new AnalysisProtocolCollection();
    	return analysisProtocolCollection;
    }

    private DataCollection generateDataCollection()
    {
    	DataCollection dataCollection = new DataCollection();
    	return dataCollection;
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
    
	public static void main(String argv[]) throws Exception
	{
		MZIdentMLGen gen = new MZIdentMLGen();
		gen.setPersonalInfo("Sangtae", "Kim", "PNNL");
        gen.writeResults(System.out);
	}
}
