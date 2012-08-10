package edu.ucsd.msjava.psi;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.Vector;
import java.util.Map.Entry;

import edu.ucsd.msjava.msdbsearch.CompactSuffixArray;
import edu.ucsd.msjava.msdbsearch.DatabaseMatch;
import edu.ucsd.msjava.msdbsearch.ScoredSpectraMap;
import edu.ucsd.msjava.msdbsearch.SearchParams;
import edu.ucsd.msjava.msgf.MSGFDBResultGenerator;
import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msscorer.SimpleDBSearchScorer;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.ModifiedAminoAcid;
import edu.ucsd.msjava.msutil.Pair;
import edu.ucsd.msjava.msutil.SpecKey;

import uk.ac.ebi.jmzidml.model.mzidml.AbstractContact;
import uk.ac.ebi.jmzidml.model.mzidml.Affiliation;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisData;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.ContactRole;
import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.DataCollection;
import uk.ac.ebi.jmzidml.model.mzidml.Enzyme;
import uk.ac.ebi.jmzidml.model.mzidml.Enzymes;
import uk.ac.ebi.jmzidml.model.mzidml.FileFormat;
import uk.ac.ebi.jmzidml.model.mzidml.FragmentationTable;
import uk.ac.ebi.jmzidml.model.mzidml.InputSpectra;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.Measure;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.jmzidml.model.mzidml.ModificationParams;
import uk.ac.ebi.jmzidml.model.mzidml.Organization;
import uk.ac.ebi.jmzidml.model.mzidml.Param;
import uk.ac.ebi.jmzidml.model.mzidml.ParamList;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.Person;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;
import uk.ac.ebi.jmzidml.model.mzidml.Role;
import uk.ac.ebi.jmzidml.model.mzidml.SearchDatabase;
import uk.ac.ebi.jmzidml.model.mzidml.SearchDatabaseRef;
import uk.ac.ebi.jmzidml.model.mzidml.SearchModification;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SourceFile;
import uk.ac.ebi.jmzidml.model.mzidml.SpectraData;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentification;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationProtocol;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.model.mzidml.Tolerance;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;

public class MZIdentMLGen {
	private MzIdentMLMarshaller m;
	
	private Person docOwner;
	private Organization org;
	private SearchDatabase searchDB;
	private ScoredSpectraMap specScanner;
	
	private final SearchParams params;
	private AminoAcidSet aaSet;
	private CompactSuffixArray sa;
	 
	private Map<edu.ucsd.msjava.msutil.Modification, Modification> modMap;
	
	private float eValueThreshold = Float.MAX_VALUE;

	// highest level objects
	private CvList cvList;
	
	private AnalysisSoftwareList analysisSoftwareList;

	// skip Provider and AuditCollection
//	private Provider provider;
//	private AuditCollection auditCollection;
	
	private SequenceCollection sequenceCollection;
	private List<DBSequence> dbSequenceList;	// list of proteins
	private List<Peptide> peptideList;			// list of peptides
	private List<PeptideEvidence> peptideEvidenceList;	// list of peptide to protein matches
	
	private AnalysisCollection analysisCollection;
	
	private AnalysisProtocolCollection analysisProtocolCollection;
	
	private DataCollection dataCollection;
	private SpectraData spectraData;
	private SpectrumIdentificationList siList;	// set of PSMs
	
	private Map<String, Peptide> pepMap;
	private Map<String, DBSequence> dbSeqMap;
	
	
	public MZIdentMLGen(SearchParams params, AminoAcidSet aaSet, CompactSuffixArray sa, ScoredSpectraMap specScanner)
	{
		m = new MzIdentMLMarshaller();
		this.params = params;
		this.aaSet = aaSet;
		this.sa = sa;
		this.specScanner = specScanner;

		modMap = new HashMap<edu.ucsd.msjava.msutil.Modification, Modification>();
		// TODO set-up modMap, SearchModification
		pepMap = new LinkedHashMap<String, Peptide>();
		dbSeqMap = new LinkedHashMap<String, DBSequence>();

		init();
	}

	private void init()
	{
		generateCvList();
		generateAnalysisSoftwareList();
		// skip Provider, AuditCollection
		initSequenceCollection();
		generateAnalysisCollection();
		generateAnalysisProtocolCollection();
		initDataCollection();
	}

	public void writeResults(PrintStream out)
	{
		out.println(m.createXmlHeader());
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
	}	
	
	public MZIdentMLGen setEValueThreshold(float eValueThreshold)
	{
		this.eValueThreshold = eValueThreshold;
		return this;
	}
	
	public Peptide addPeptide(int index, int length, String pepStr)
	{
//		char pre = sa.getSequence().getCharAt(index);
//		char post = sa.getSequence().getCharAt(index+length-1);
//		String pvID = pre+pepStr+post;
		
		Peptide mzidPeptide = pepMap.get(pepStr);
		
		if(mzidPeptide != null)
			return mzidPeptide;
				
		// new peptide variant
		mzidPeptide = new Peptide();
		List<Modification> modList = mzidPeptide.getModification();
		edu.ucsd.msjava.msutil.Peptide peptide = aaSet.getPeptide(pepStr); 
		StringBuffer unmodPepStr = new StringBuffer();
		int location = 0;
		for(edu.ucsd.msjava.msutil.AminoAcid aa : peptide)
		{
			unmodPepStr.append(aa.getUnmodResidue());
			if(aa.isModified())
			{
				Modification mod = new Modification();
				mod.setLocation(location);
				ModifiedAminoAcid modAA = (ModifiedAminoAcid)aa;
				mod.setMonoisotopicMassDelta(modAA.getModification().getAccurateMass());
				modList.add(mod);
			}
			location++;
		}

		mzidPeptide.setId(pepStr);
		pepMap.put(pepStr, mzidPeptide);
		
		return mzidPeptide;
	}
	
	public synchronized void addSpectrumIdentificationResults(Map<SpecKey,PriorityQueue<DatabaseMatch>> specKeyDBMatchMap)
	{
		Map<Integer,PriorityQueue<DatabaseMatch>> specIndexDBMatchMap = new HashMap<Integer,PriorityQueue<DatabaseMatch>>();
		
		Iterator<Entry<SpecKey, PriorityQueue<DatabaseMatch>>> itr = specKeyDBMatchMap.entrySet().iterator();
		int numPeptidesPerSpec = params.getNumMatchesPerSpec();
		while(itr.hasNext())
		{
			Entry<SpecKey, PriorityQueue<DatabaseMatch>> entry = itr.next();
			SpecKey specKey = entry.getKey();
			PriorityQueue<DatabaseMatch> matchQueue = entry.getValue();
			if(matchQueue == null || matchQueue.size() == 0)
				continue;
			
			int specIndex = specKey.getSpecIndex();
			PriorityQueue<DatabaseMatch> existingQueue = specIndexDBMatchMap.get(specIndex);
			if(existingQueue == null)
			{
				existingQueue = new PriorityQueue<DatabaseMatch>(numPeptidesPerSpec, new DatabaseMatch.SpecProbComparator());
				specIndexDBMatchMap.put(specIndex, existingQueue);
			}
			
			for(DatabaseMatch match : matchQueue)
			{
				if(existingQueue.size() < numPeptidesPerSpec)
				{
					existingQueue.add(match);
				}
				else if(existingQueue.size() >= numPeptidesPerSpec)
				{
					if(match.getSpecEValue() < existingQueue.peek().getSpecEValue())
					{
						existingQueue.poll();
						existingQueue.add(match);
					}
				}
			}
		}		
		
		Iterator<Entry<Integer, PriorityQueue<DatabaseMatch>>> itr2 = specIndexDBMatchMap.entrySet().iterator();
		while(itr2.hasNext())
		{
			Entry<Integer, PriorityQueue<DatabaseMatch>> entry = itr2.next();
			int specIndex = entry.getKey();
			
			PriorityQueue<DatabaseMatch> matchQueue = entry.getValue();
			if(matchQueue == null || matchQueue.size() == 0)
				continue;

			
			Pair<String,Float> idAndSpecIndexPair = specScanner.getIDAndPrecursorMzPair(specIndex);
			String id = idAndSpecIndexPair.getFirst();
			float precursorMz = idAndSpecIndexPair.getSecond();

			SpectrumIdentificationResult sir = new SpectrumIdentificationResult();
			sir.setId(id);
			sir.setSpectraData(spectraData);
			
			ArrayList<DatabaseMatch> matchList = new ArrayList<DatabaseMatch>(matchQueue);
			int rank = 0;
			
			for(int i=matchList.size()-1; i>=0; --i)
			{
				++rank;
				DatabaseMatch match = matchList.get(i);
				
				if(match.getDeNovoScore() < 0)
					break;
				
				int pepIndex = match.getIndex();	// Position of preAA
				int length = match.getLength();		// Peptide length + 2
				int charge = match.getCharge();
				String peptideStr = match.getPepSeq();

				float peptideMass = match.getPeptideMass();
				float theoMass = peptideMass + (float)Composition.H2O;
				float theoMz = theoMass/charge + (float)Composition.H;
				
				int score = match.getScore();
				double specEValue = match.getSpecEValue();
				int numPeptides = sa.getNumDistinctPeptides(params.getEnzyme() == null ? length-2 : length-1);
				double eValue = specEValue*numPeptides;
				
				String specEValueStr;
				if(specEValue < Float.MIN_NORMAL)
					specEValueStr = String.valueOf(specEValue);
				else
					specEValueStr = String.valueOf((float)specEValue);
				
				String eValueStr;
				if(specEValue < Float.MIN_NORMAL)
					eValueStr = String.valueOf(eValue);
				else
					eValueStr = String.valueOf((float)eValue);

				SpectrumIdentificationItem sii = new SpectrumIdentificationItem();

				sii.setChargeState(charge);
				sii.setExperimentalMassToCharge(precursorMz);
				sii.setCalculatedMassToCharge((double) theoMz);
				sii.setPeptide(addPeptide(pepIndex, length, peptideStr));
				sii.setRank(rank);
				sii.setPassThreshold(eValue <= eValueThreshold);
				sii.setId(Constants.siiID+specIndex+"_"+rank);
					
				List<CvParam> cvList = sii.getCvParam();
				List<UserParam> userList = sii.getUserParam();

				CvParam rawScoreCV = Constants.makeCvParam("MS:1002049", "MS-GF:RawScore");
				rawScoreCV.setValue(String.valueOf(score));
				cvList.add(rawScoreCV);

				CvParam deNovoScoreCV = Constants.makeCvParam("MS:1002050", "MS-GF:DeNovoScore");
				deNovoScoreCV.setValue(String.valueOf(match.getDeNovoScore()));
				cvList.add(deNovoScoreCV);

				CvParam specEValueCV = Constants.makeCvParam("MS:1002052", "MS-GF:SpecEValue");
				specEValueCV.setValue(specEValueStr);
				cvList.add(specEValueCV);

				CvParam eValueCV = Constants.makeCvParam("MS:1002052", "MS-GF:EValue");
				eValueCV.setValue(eValueStr);
				cvList.add(eValueCV);
				
				UserParam isotopeErrorParam = Constants.makeUserParam("Isotope error");
				float expMass = (precursorMz - (float)Composition.H)*charge;
				int isotopeError = NominalMass.toNominalMass(expMass) - NominalMass.toNominalMass(theoMass);
				isotopeErrorParam.setValue(String.valueOf(isotopeError));
				userList.add(isotopeErrorParam);
				
				sir.getSpectrumIdentificationItem().add(sii);
			}
			siList.getSpectrumIdentificationResult().add(sir);
		}
	}	
	
	private void generateCvList()
	{
		cvList = new CvList();
		List<Cv> localCvList = cvList.getCv();
		localCvList.add(Constants.psiCV);
		localCvList.add(Constants.unimodCV);
		localCvList.add(Constants.unitCV);
	}

	private void generateAnalysisSoftwareList()
	{
		analysisSoftwareList = new AnalysisSoftwareList();
		List<AnalysisSoftware> analysisSoftwares = analysisSoftwareList.getAnalysisSoftware();
		analysisSoftwares.add(Constants.msgfPlus);
	}
	
	private void initSequenceCollection()
	{
		sequenceCollection = new SequenceCollection();
		dbSequenceList = sequenceCollection.getDBSequence();
		peptideList = sequenceCollection.getPeptide();
		peptideEvidenceList = sequenceCollection.getPeptideEvidence();
	}
	
	private void generateAnalysisCollection()
	{
		analysisCollection = new AnalysisCollection();
//		List<SpectrumIdentification> specIdentList = analysisCollection.getSpectrumIdentification();
//		SpectrumIdentification specIdent = new SpectrumIdentification();
//		specIdent.setId(Constants.specIdentID);
//
//		List<SpectrumIdentification> siList = analysisCollection.getSpectrumIdentification();
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
//		List<InputSpectra> inputSpecList = specIdent.getInputSpectra();
//		InputSpectra inputSpec = new InputSpectra();
//		inputSpec.setSpectraData(spectraData);
//		inputSpecList.add(inputSpec);
//
//		specIdentList.add(specIdent);
	}	
	
	private void generateAnalysisProtocolCollection()
	{
		AnalysisProtocolCollectionGen apcGen = new AnalysisProtocolCollectionGen(params, aaSet);
		analysisProtocolCollection = apcGen.getAnalysisProtocolCollection();
	}
	
	private void initDataCollection()
	{
		dataCollection = new DataCollection();

		// Inputs
		Inputs inputs = new Inputs();
		// source file
		// search database
		// spectra data
		spectraData = new SpectraData();
		spectraData.setId(Constants.spectraDataID);
		spectraData.setName("SpectrumName");	// TODO: set up name
		// TODO: set up file format
		inputs.getSpectraData().add(spectraData);
		dataCollection.setInputs(inputs);
		
		// AnalysisData
		AnalysisData analysisData = new AnalysisData();
		dataCollection.setAnalysisData(analysisData);
		
		siList = new SpectrumIdentificationList();
		siList.setId(Constants.siiListID);
		
		analysisData.getSpectrumIdentificationList().add(siList);
	}
	

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
//
//	private void generateProvider()
//	{
//		if(docOwner == null)
//			return null;
//
//		provider = new Provider();
//		provider.setId(Constants.providerID);
//
//		ContactRole contactRole = new ContactRole();
//		contactRole.setContact(docOwner);

//		Role role = new Role();
//		role.setCvParam(makeCvParam("MS:1001271","researcher",psiCV));
//		contactRole.setRole(role);
//		provider.setContactRole(contactRole);
//
//		return provider;
//	}
//
//
//	private void generateAuditCollection()
//	{
//		if(docOwner == null || org == null)
//			return null;
//
//		auditCollection = new AuditCollection();
//		List<AbstractContact> contactList = auditCollection.getPersonOrOrganization();
//
//		contactList.add(docOwner);
//		contactList.add(org);
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
}
