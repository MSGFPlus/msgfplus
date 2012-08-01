//package psi;
//
//import java.io.PrintStream;
//import java.util.List;
//import java.util.Vector;
//
//import uk.ac.ebi.jmzidml.model.mzidml.AbstractContact;
//import uk.ac.ebi.jmzidml.model.mzidml.Affiliation;
//import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
//import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
//import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
//import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
//import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
//import uk.ac.ebi.jmzidml.model.mzidml.ContactRole;
//import uk.ac.ebi.jmzidml.model.mzidml.Cv;
//import uk.ac.ebi.jmzidml.model.mzidml.CvList;
//import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
//import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
//import uk.ac.ebi.jmzidml.model.mzidml.DataCollection;
//import uk.ac.ebi.jmzidml.model.mzidml.Enzyme;
//import uk.ac.ebi.jmzidml.model.mzidml.Enzymes;
//import uk.ac.ebi.jmzidml.model.mzidml.FileFormat;
//import uk.ac.ebi.jmzidml.model.mzidml.FragmentationTable;
//import uk.ac.ebi.jmzidml.model.mzidml.InputSpectra;
//import uk.ac.ebi.jmzidml.model.mzidml.Measure;
//import uk.ac.ebi.jmzidml.model.mzidml.Modification;
//import uk.ac.ebi.jmzidml.model.mzidml.ModificationParams;
//import uk.ac.ebi.jmzidml.model.mzidml.Organization;
//import uk.ac.ebi.jmzidml.model.mzidml.Param;
//import uk.ac.ebi.jmzidml.model.mzidml.ParamList;
//import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
//import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
//import uk.ac.ebi.jmzidml.model.mzidml.Person;
//import uk.ac.ebi.jmzidml.model.mzidml.Provider;
//import uk.ac.ebi.jmzidml.model.mzidml.Role;
//import uk.ac.ebi.jmzidml.model.mzidml.SearchDatabase;
//import uk.ac.ebi.jmzidml.model.mzidml.SearchDatabaseRef;
//import uk.ac.ebi.jmzidml.model.mzidml.SearchModification;
//import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
//import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentification;
//import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
//import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationProtocol;
//import uk.ac.ebi.jmzidml.model.mzidml.Tolerance;
//import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
//import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
//
//public class MZIdentMLGen {
//	private MzIdentMLMarshaller m;
//	private Person docOwner;
//	private Organization org;
//	private SearchDatabase searchDB;
//
//	public MZIdentMLGen()
//	{
//		m = new MzIdentMLMarshaller();
//		initCVs();
//	}
//
//	public void writeResults(PrintStream out)
//	{
//		out.println(m.createXmlHeader());
//		out.println(m.createMzIdentMLStartTag("12345"));
//		m.marshal(generateCvList(), out);
//		out.println();
//		m.marshal(generateAnalysisSoftwareList(), out);
//		out.println();
//
//		Provider provider = generateProvider();
//		if(provider != null)
//		{
//			m.marshal(generateProvider(), out);
//			out.println();
//		}
//		AuditCollection auditCollection = generateAuditCollection();
//		if(auditCollection != null)
//		{
//			m.marshal(generateAuditCollection(), out);
//			out.println();
//		}
//
//		m.marshal(generateSequenceCollection(), out);
//		out.println();
//
//		m.marshal(generateAnalysisCollection(), out);
//		out.println();
//
//		m.marshal(generateAnalysisProtocolCollection(), out);
//		out.println();
//
//		m.marshal(generateDataCollection(), out);
//		out.println();
//	}
//
//	public void setDatabses()
//	{
//		searchDB = new SearchDatabase();
//		searchDB.setId(searchDBID);
//		searchDB.setNumDatabaseSequences(100L);
//
//		UserParam param = new UserParam();
//		param.setName("DatabaseName");
//		Param tempParam = new Param();
//		tempParam.setParam(param);
//		searchDB.setDatabaseName(tempParam);
//
//		searchDB.setLocation("DBLocation");
//		FileFormat ff = new FileFormat();
//		ff.setCvParam(makeCvParam("MS:1001348","FASTA format",psiCV));
//		searchDB.setFileFormat(ff);   
//	}
//
//	public MZIdentMLGen setPersonalInfo(
//			String firstName,
//			String lastName,
//			String affiliationName
//			)
//	{
//		docOwner = new Person();
//		docOwner.setId("PERSON_DOC_OWNER");
//		docOwner.setFirstName(firstName);
//		docOwner.setLastName(lastName);
//
//		org = new Organization();
//		org.setId("ORG_DOC_OWNER");
//		org.setName(affiliationName);
//
//		Affiliation aff = new Affiliation();
//		aff.setOrganization(org);
//		docOwner.getAffiliation().add(aff);
//
//		return this;
//	}
//
//	private CvList generateCvList()
//	{
//		CvList cvList = new CvList();
//		List<Cv> localCvList = cvList.getCv();
//
//		Cv unimodCV = new Cv();
//		unimodCV.setUri("http://www.unimod.org/obo/unimod.obo");
//		unimodCV.setId(unimodID);
//		unimodCV.setFullName("UNIMOD");
//
//		Cv unitCV = new Cv();
//		unitCV.setUri("http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo");
//		unitCV.setId(unitCvID);
//		unitCV.setFullName("UNIT-ONTOLOGY");
//
//		localCvList.add(psiCV);
//		localCvList.add(unimodCV);
//		localCvList.add(unitCV);
//
//		return cvList;
//	}
//
//	private AnalysisSoftwareList generateAnalysisSoftwareList()
//	{
//		AnalysisSoftwareList analysisSoftwareList = new AnalysisSoftwareList();
//		List<AnalysisSoftware> analysisSoftwares = analysisSoftwareList.getAnalysisSoftware();
//		AnalysisSoftware analysisSoftware = new AnalysisSoftware();
//		analysisSoftware.setName("MS-GF+");
//		Param tempParam = new Param();
//		tempParam.setParam(makeCvParam("MS:XXXXXXX","MS-GF+",psiCV));
//		analysisSoftware.setSoftwareName(tempParam);
//		analysisSoftware.setId(analysisSoftID);
//		analysisSoftware.setVersion("1.0");
//		analysisSoftwares.add(analysisSoftware);
//
//		return analysisSoftwareList;    	
//	}
//
//	// TODO: check whether Provider information can be extracted from mzML
//	private Provider generateProvider()
//	{
//		if(docOwner == null)
//			return null;
//
//		Provider provider = new Provider();
//		provider.setId("PROVIDER");
//
//		ContactRole contactRole = new ContactRole();
//		contactRole.setContact(docOwner);
//
//		Role role = new Role();
//		role.setCvParam(makeCvParam("MS:1001271","researcher",psiCV));
//		contactRole.setRole(role);
//		provider.setContactRole(contactRole);
//
//		return provider;
//	}
//
//
//	private AuditCollection generateAuditCollection()
//	{
//		if(docOwner == null || org == null)
//			return null;
//
//		AuditCollection auditCollection = new AuditCollection();
//		List<AbstractContact> contactList = auditCollection.getPersonOrOrganization();
//
//		contactList.add(docOwner);
//		contactList.add(org);
//
//		return auditCollection;
//	}
//
//	private SequenceCollection generateSequenceCollection()
//	{
//		SequenceCollection sequenceCollection = new SequenceCollection();
//		// Set List<DBSequence>, List<Peptide>, List<PeptideEvidence>
//
//		String pepSeq = "PEPTIDE";
//		String protAccession = "PROTEINAcc";
//		String dbName = "DBName.fasta";
//
//		//  Set List<DBSequence>
//		List<DBSequence> dbSequenceList = sequenceCollection.getDBSequence();
//		DBSequence dbSeq = new DBSequence();
//		dbSeq.setAccession(protAccession);
//		dbSeq.setLength(100);
//		dbSeq.setId("DBSeq_" + protAccession);
//		dbSeq.setSearchDatabase(searchDB);
//		List<CvParam> dbCvParamList = dbSeq.getCvParam();
//		CvParam protDescCV = makeCvParam("MS:1001088","protein description", psiCV);
//		protDescCV.setName("protein description");
//		dbCvParamList.add(protDescCV);
//		dbSequenceList.add(dbSeq);
//
//		// Set List<Peptide>
//		List<Peptide> peptideList = sequenceCollection.getPeptide();
//
//		// set peptide
//		Peptide mzidPep = new Peptide();
//		mzidPep.setPeptideSequence(pepSeq);
//		mzidPep.setId("PepID");
//		List<Modification> allMods = mzidPep.getModification();
//
//		// set modification
//		Modification mzidmod = new Modification();
//		double monoMass = msutil.Modification.get("Oxidation").getAccurateMass();
//
//		// set modification CV params
//		List<CvParam> paramList = mzidmod.getCvParam();
//		CvParam modParam = new CvParam();
//		int unimodRecordID = 35;
//		modParam.setAccession("UNIMOD:" + unimodRecordID);	// TODO: use Unimod
//		modParam.setCv(unimodCV);
//		modParam.setName("Oxidation");	// TODO: use unimod name
//		paramList.add(modParam);
//
//		mzidmod.setMonoisotopicMassDelta(monoMass);
//		mzidmod.setLocation(1);   // 1-based
//
//		// add modification to the peptide
//		allMods.add(mzidmod);
//
//		// add mzidPep to peptideList
//		peptideList.add(mzidPep);
//
//		// Set List<PeptideEvidence>
//		List<PeptideEvidence> peptideEvidenceList = sequenceCollection.getPeptideEvidence();
//
//		PeptideEvidence pepEvid = new PeptideEvidence();
//		pepEvid.setStart(0);
//		pepEvid.setEnd(15);
//		pepEvid.setPre("K");
//		pepEvid.setPost("G");
//		pepEvid.setId("PepEvID");
//		pepEvid.setIsDecoy(false);
//		pepEvid.setPeptide(mzidPep);
//		pepEvid.setDBSequence(dbSeq);
//
//		peptideEvidenceList.add(pepEvid);
//
//		return sequenceCollection;
//	}
//
//	private AnalysisCollection generateAnalysisCollection()
//	{
//		AnalysisCollection analysisCollection = new AnalysisCollection();
//		List<SpectrumIdentification> specIdentList = analysisCollection.getSpectrumIdentification();
//		SpectrumIdentification specIdent = new SpectrumIdentification();
//		specIdent.setId("SpecIdentID");
//
//		SpectrumIdentificationList siList = new SpectrumIdentificationList();
//		siList.setId(siiListID);
//
//		
//		FragmentationTable fragTable = new FragmentationTable();
//		List<Measure> measureList = fragTable.getMeasure();
//		Measure mzMeasure = new Measure();
//		mzMeasure.setId(measureMzID);
//		List<CvParam> cvParamList = mzMeasure.getCvParam();
//		cvParamList.add(makeCvParam("MS:1001225","product ion m/z", psiCV, "MS:1000040", "m/z", psiCV));
//		Measure intMeasure = new Measure();
//		intMeasure.setId(measureIntID);
//		cvParamList = intMeasure.getCvParam();
//		cvParamList.add(makeCvParam("MS:1001226","product ion intensity",psiCV,"MS:1000131", "number of counts", psiCV));
//		Measure errorMeasure = new Measure();
//		errorMeasure.setId(measureErrorID);
//		cvParamList = errorMeasure.getCvParam();
//		cvParamList.add(makeCvParam("MS:1001227","product ion m/z error",psiCV,"MS:1000040","m/z", psiCV));
//
//		measureList.add(mzMeasure);
//		measureList.add(intMeasure);
//		measureList.add(errorMeasure);
//
//		siList.setFragmentationTable(fragTable);
//
//		specIdent.setSpectrumIdentificationList(siList);
//		specIdent.setSpectrumIdentificationProtocol(siProtocol);
//		List<SearchDatabaseRef> searchDBRefList = specIdent.getSearchDatabaseRef();
//		SearchDatabaseRef searchDBRef = new SearchDatabaseRef();
//		searchDBRef.setSearchDatabase(searchDB);
//		searchDBRefList.add(searchDBRef);
//
//
//		List<InputSpectra> inputSpecList = specIdent.getInputSpectra();
//		InputSpectra inputSpec = new InputSpectra();
//		inputSpec.setSpectraData(spectraData);
//		inputSpecList.add(inputSpec);
//
//		specIdentList.add(specIdent);
//
//		return analysisCollection;
//	}
//	
//	private void generateSearchProtocol()
//	{
//		SpectrumIdentificationProtocol siProtocol = new SpectrumIdentificationProtocol();
//        siProtocol.setId(siProtocolID);
//        siProtocol.setAnalysisSoftware(analysisSoftware);
//
//        //<cvParam accession="MS:1001083" name="ms-ms search" cvRef="PSI-MS"/>
//        Param tempParam = new Param();
//        tempParam.setParam(makeCvParam("MS:1001083","ms-ms search",psiCV));
//        siProtocol.setSearchType(tempParam);
//
//
//        //List<CvParam> cvParamList = siProtocol.getAdditionalSearchCvParams();
//        ParamList paramList = siProtocol.getAdditionalSearchParams();
//        if(paramList == null ){
//            paramList = new ParamList();
//            siProtocol.setAdditionalSearchParams(paramList);
//        }
//        List<CvParam> cvParamList = paramList.getCvParam();
//
//        int msSearchType = settings.MSSearchSettings_precursorsearchtype.MSSearchType;    //with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact
//
//        if(msSearchType==0){
//            cvParamList.add(makeCvParam("MS:1001211","parent mass type mono",psiCV));
//        }
//        else if(msSearchType==1){
//            cvParamList.add(makeCvParam("MS:1001212","parent mass type average",psiCV));
//        }
//        else if(msSearchType==2){
//            cvParamList.add(makeCvParam("MS:No_acc","monoisotopic N15",psiCV));
//            System.out.println("Warning: No CV term for monoisotopic N15 search type");
//        }
//        else if(msSearchType==3){
//            cvParamList.add(makeCvParam("MS:No_acc","exact",psiCV));
//            System.out.println("Warning: No CV term for exact mass search type");
//        }
//        else{
//            System.out.println("Error search type not recognised");
//
//        }
//
//        int prodSearchType = settings.MSSearchSettings_productsearchtype.MSSearchType;    //with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact
//
//        if(prodSearchType==0){
//            cvParamList.add(makeCvParam("MS:1001256","fragment mass type mono",psiCV));
//        }
//        else if(prodSearchType==1){
//            cvParamList.add(makeCvParam("MS:1001255","fragment mass type average",psiCV));
//        }
//        else if(prodSearchType==2){
//            cvParamList.add(makeCvParam("MS:No_acc","monoisotopic N15",psiCV));
//            System.out.println("Warning: No CV term for monoisotopic N15 search type");
//        }
//        else if(prodSearchType==3){
//            cvParamList.add(makeCvParam("MS:No_acc","exact",psiCV));
//            System.out.println("Warning: No CV term for exact mass search type");
//        }
//        else{
//            System.out.println("Error search type not recognised");
//        }
//
//
//        ModificationParams modParams = new ModificationParams();
//        List<SearchModification> searchModList = modParams.getSearchModification();
//
//        for(int fixedMod : settings.MSSearchSettings_fixed.MSMod){
//            OmssaModification omod = intToModMap.get(fixedMod);
//            double monoMass = omod.getModMonoMass();
//            Vector<String> residues = omod.getModResidues();
//            Boolean isMono = true;
//            ModT unimod = unimodDoc.getModByMass(monoMass, unimodMassError, isMono, residues);
//            //getModByMass(Double testMass, Double massError, Boolean isMono, Vector<String> residues)
//
//            SearchModification searchMod = new SearchModification();
//
//            //One example for variable Oxidation on M
//            searchMod.setFixedMod(true);
//            //ModParam modParam = new ModParam();
//            //modParam.setMassDelta(new Float(monoMass));
//            //modParam.setCvParam();
//            //searchMod.setModParam(modParam);
//            List<CvParam> modCvParamList = searchMod.getCvParam();
//            modCvParamList.add(makeCvParam("UNIMOD:" + unimod.getRecordId(),unimod.getTitle(),unimodCV));
//            searchMod.setMassDelta(new Float(monoMass));
//            List<String> residueList = searchMod.getResidues();
//
//            for(String residue : residues){
//                residueList.add(residue);
//            }
//
//            searchModList.add(searchMod);
//        }
//
//        for(int varMod : settings.MSSearchSettings_variable.MSMod){
//            OmssaModification omod = intToModMap.get(varMod);
//            double monoMass = omod.getModMonoMass();
//            Vector<String> residues = omod.getModResidues();
//            Boolean isMono = true;
//            ModT unimod = unimodDoc.getModByMass(monoMass, unimodMassError, isMono, residues);
//
//            SearchModification searchMod = new SearchModification();
//
//            //One example for variable Oxidation on M
//            searchMod.setFixedMod(false);
//            //ModParam modParam = new ModParam();
//            //modParam.setMassDelta(new Float(monoMass));
//            //modParam.setCvParam(makeCvParam("UNIMOD:" + unimod.getRecordId(),unimod.getTitle(),unimodCV));
//            //searchMod.setModParam(modParam);
//            searchModList.add(searchMod);
//
//
//            List<CvParam> modCvParamList = searchMod.getCvParam();
//            modCvParamList.add(makeCvParam("UNIMOD:" + unimod.getRecordId(),unimod.getTitle(),unimodCV));
//            searchMod.setMassDelta(new Float(monoMass));
//            List<String> residueList = searchMod.getResidues();
//
//            for(String residue : residues){
//                residueList.add(residue);
//            }
//        }
//
//        siProtocol.setModificationParams(modParams);
//        
//
//                /*
//            <ModificationParams>
//                <SearchModification fixedMod="false">
//                          <ModParam massDelta="15.994919" residues="M">
//                            <cvParam accession="UNIMOD:35" name="Oxidation" cvRef="UNIMOD" />
//                         </ModParam>
//                </SearchModification>
//                */
//
//        /*
//            <Enzymes independent="0">
//                <Enzyme id="ENZ_1" CTermGain="OH" NTermGain="H" missedCleavages="1" semiSpecific="0">
//                <EnzymeName>
//                <cvParam accession="MS:1001251" name="Trypsin" cvRef="PSI-MS" />
//                </EnzymeName>
//                </Enzyme>
//             </Enzymes>
//*/
//
//
//        Enzymes enzymes = siProtocol.getEnzymes();
//
//        if(enzymes==null){
//            enzymes = new Enzymes();
//            siProtocol.setEnzymes(enzymes);
//        }
//        enzymes.setIndependent(false);
//
//        List<Enzyme> enzymeList = enzymes.getEnzyme();
//
//        List<Integer> msEnzymeList = settings.MSSearchSettings_enzyme.MSEnzymes;
//
//        for (Integer msEnzyme : msEnzymeList){
//            OmssaEnumerators omssaEnums = new OmssaEnumerators();
//            Enzyme enzyme = getEnzyme(omssaEnums.getEnzymeAsText(msEnzyme),settings.MSSearchSettings_missedcleave);
//            enzymeList.add(enzyme);
//        }
//
//
//        Tolerance fragTol = new Tolerance();
//        Tolerance parTol = new Tolerance();
//
//        List<CvParam> fragCvList = fragTol.getCvParam();
//        CvParam fragCvPlus = getCvParamWithMassUnits(true);
//        CvParam fragCvMinus = getCvParamWithMassUnits(true);
//
//
//        /*
//         <FragmentTolerance>
//                            <cvParam accession="MS:1001412" name="search tolerance plus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
//                            <cvParam accession="MS:1001413" name="search tolerance minus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
//                    </FragmentTolerance>
//                    <ParentTolerance>
//                            <cvParam accession="MS:1001412" name="search tolerance plus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
//                            <cvParam accession="MS:1001413" name="search tolerance minus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
//                    </ParentTolerance>
//          */
//
//        fragCvPlus.setAccession("MS:1001412");
//        fragCvPlus.setName("search tolerance plus value");
//        fragCvMinus.setAccession("MS:1001413");
//        fragCvMinus.setName("search tolerance minus value");
//        fragCvPlus.setValue(""+settings.MSSearchSettings_msmstol);    
//        fragCvMinus.setValue(""+settings.MSSearchSettings_msmstol);
//        fragCvList.add(fragCvPlus);
//        fragCvList.add(fragCvMinus);
//
//        List<CvParam> parCvList = parTol.getCvParam();
//        CvParam parCvPlus = getCvParamWithMassUnits(true);
//        CvParam parCvMinus = getCvParamWithMassUnits(true);
//
//        parCvPlus.setAccession("MS:1001412");
//        parCvPlus.setName("search tolerance plus value");
//        parCvMinus.setAccession("MS:1001413");
//        parCvMinus.setName("search tolerance minus value");
//        parCvPlus.setValue(""+settings.MSSearchSettings_peptol); 
//        parCvMinus.setValue(""+settings.MSSearchSettings_peptol);
//        parCvList.add(parCvPlus);
//        parCvList.add(parCvMinus);
//
//        siProtocol.setFragmentTolerance(fragTol);
//        siProtocol.setParentTolerance(parTol);
//
//        // siProtocol.getThresholdCvParams();
//        ParamList sip_paramList = siProtocol.getThreshold();
//        if(sip_paramList == null ){
//            sip_paramList = new ParamList();
//            siProtocol.setThreshold(sip_paramList);
//        }
//        cvParamList =sip_paramList.getCvParam();
//
//        cvParamList.add(makeCvParam("MS:1001494","no threshold",psiCV));
//        //<cvParam accession="MS:1001494" name="no threshold" cvRef="PSI-MS" />
//        
//        sipList.add(siProtocol);		
//	}
//
//	private AnalysisProtocolCollection generateAnalysisProtocolCollection()
//	{
//		AnalysisProtocolCollection analysisProtocolCollection = new AnalysisProtocolCollection();
//		return analysisProtocolCollection;
//	}
//
//	private DataCollection generateDataCollection()
//	{
//		DataCollection dataCollection = new DataCollection();
//		return dataCollection;
//	}
//
//	public static void main(String argv[]) throws Exception
//	{
//		MZIdentMLGen gen = new MZIdentMLGen();
//		gen.setPersonalInfo("Sangtae", "Kim", "PNNL");
//		gen.writeResults(System.out);
//	}
//}
