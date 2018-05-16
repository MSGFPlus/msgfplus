package edu.ucsd.msjava.mzid;

import edu.ucsd.msjava.msdbsearch.*;
import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.msutil.Modification.Location;
import edu.ucsd.msjava.ui.MSGFPlus;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.model.mzidml.FileFormat;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;

import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.stream.StreamSource;
import java.io.*;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.*;

public class MZIdentMLGen {
    private MzIdentMLMarshaller m;

//	private Person docOwner;
//	private Organization org;

    private final SearchParams params;
    private AminoAcidSet aaSet;
    private CompactSuffixArray sa;
    private SpectraAccessor specAcc;
    private final int ioIndex;

    private float eValueThreshold = Float.MAX_VALUE;

    // highest level objects
    private CvList cvList;

    private AnalysisSoftwareList analysisSoftwareList;

    // skip Provider and AuditCollection
//	private Provider provider;
//	private AuditCollection auditCollection;

    private SequenceCollection sequenceCollection;
    private List<DBSequence> dbSequenceList;    // list of proteins
    private List<Peptide> peptideList;            // list of peptides
    private List<PeptideEvidence> peptideEvidenceList;    // list of peptide to protein matches

    private AnalysisCollection analysisCollection;

    private AnalysisProtocolCollection analysisProtocolCollection;
    private SpectrumIdentificationProtocol siProtocol;

    private DataCollection dataCollection;
    private SpectraData spectraData;
    private SpectrumIdentificationList siList;    // set of PSMs
    private SearchDatabase searchDatabase;

    private Map<Integer, DBSequence> dbSeqMap;
    private Map<Integer, Boolean> isDecoyMap;
    private Map<String, Peptide> pepMap;
    private Map<String, PeptideEvidence> pepEvMap;
    private Map<String, List<PeptideEvidenceRef>> evRefListMap;

    private AnalysisProtocolCollectionGen apcGen;

    public static final String encoding = "UTF-8";
    public static final Charset charset = StandardCharsets.UTF_8;

    public MZIdentMLGen(SearchParams params, AminoAcidSet aaSet, CompactSuffixArray sa, SpectraAccessor specAcc, int ioIndex) {
        m = new MzIdentMLMarshaller();
        this.params = params;
        this.aaSet = aaSet;
        this.sa = sa;
        this.specAcc = specAcc;
        this.ioIndex = ioIndex;

        dbSeqMap = new HashMap<Integer, DBSequence>();
        isDecoyMap = new HashMap<Integer, Boolean>();
        pepMap = new LinkedHashMap<String, Peptide>();
        pepEvMap = new LinkedHashMap<String, PeptideEvidence>();
        evRefListMap = new HashMap<String, List<PeptideEvidenceRef>>();

        init();
    }

    private void init() {
        generateCvList();
        generateAnalysisSoftwareList();
        generateAnalysisProtocolCollection();
        // skip Provider, AuditCollection
        initSequenceCollection();
        initDataCollection();
        generateAnalysisCollection();
    }

    /**
     * Output the results to the specified mzIdentML file.
     *
     * @param file
     */
    public void writeResults(File file) {
        // Original function - produces e.g. <cvParam ...></cvParam> after updating jmzidentml from 1.1.3 to 1.2.0 or newer
        /*
        jmzidentml version notes: 
        1.1.3 & 1.1.4 produce good output, use self-closing tags, use 4 space indents
        1.1.5-1.1.9 produce bad output (missing elements), use self-closing tags, 2 space indents
        1.2.0+ produce good output, use closing tags instead of self-closing tags, 2 space indents
        */
        //writeResultsOriginal(file);
        //writeResultsCleanerTempString(file);
        writeResultsCleanerTempFile(file);
    }

    /**
     * Original writeResults; after a dependency update, this is no longer writing self-closing tags, resulting in larger mzid files.
     *
     * @param file
     */
    public void writeResultsOriginal(File file) {
        OutputStream os = null;
        try {
            os = new FileOutputStream(file);
        } catch (FileNotFoundException e) {
            System.out.println("Could not find file \"" + file.getAbsolutePath() + "\" to write to. Writing to console...");
            os = System.out;
        }
        OutputStreamWriter out = new OutputStreamWriter(os, charset);
        writeResults(out);
    }

    /**
     * Single function to handle the writing of the results to a Writer
     *
     * @param out
     */
    private void writeResults(Writer out) {
        try {
            //out.write(m.createXmlHeader());
            out.write("<?xml version=\"1.0\" encoding=\"" + encoding + "\"?>\n");
            out.write(m.createMzIdentMLStartTag("MS-GF+") + "\n");
            m.marshal(cvList, out);
            out.write("\n");
            m.marshal(analysisSoftwareList, out);
            out.write("\n");

            m.marshal(sequenceCollection, out);
            out.write("\n");

            m.marshal(analysisCollection, out);
            out.write("\n");

            m.marshal(analysisProtocolCollection, out);
            out.write("\n");

            m.marshal(dataCollection, out);
            out.flush();
            out.write("\n");

            out.write(m.createMzIdentMLClosingTag() + "\n");
            out.flush();
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Output results to a StringWriter, then pipe that string through a transform to make empty elements into self-closing elements
     *
     * @param file
     */
    public void writeResultsCleanerTempString(File file) {
        OutputStream os = null;
        try {
            os = new FileOutputStream(file);
        } catch (FileNotFoundException e) {
            System.out.println("Could not find file \"" + file.getAbsolutePath() + "\" to write to. Writing to console...");
            os = System.out;
        }
        try {
            StringWriter out = new StringWriter();
            writeResults(out);

            OutputStreamWriter out2 = new OutputStreamWriter(os, charset);
            StringReader sin = new StringReader(out.toString());
            TransformerFactory transFactory = TransformerFactory.newInstance();
            Transformer transformer = transFactory.newTransformer();
            transformer.transform(new StreamSource(sin), new StreamResult(out2));
            out2.flush();
            out2.close();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (TransformerException e) {
            e.printStackTrace();
        }
    }

    /**
     * Output results to a temp file, then pipe that file through a transform to make empty elements into self-closing elements
     *
     * @param file
     */
    public void writeResultsCleanerTempFile(File file) {
        FileInputStream is = null;
        FileOutputStream os = null;
        String newFileName = file.getAbsolutePath();
        int decimal = newFileName.lastIndexOf(".");
        String ext = newFileName.substring(decimal + 1);
        newFileName = newFileName.substring(0, decimal) + "_dirty." + ext;
        File tempFile = new File(newFileName);

        writeResultsOriginal(tempFile);
        try {
            is = new FileInputStream(tempFile);
            os = new FileOutputStream(file);

            InputStreamReader tempIn = new InputStreamReader(is, charset);
            OutputStreamWriter out2 = new OutputStreamWriter(os, charset);
            // Manually output the xml declaration so that it is on its own line
            out2.write("<?xml version=\"1.0\" encoding=\"" + encoding + "\"?>\n");
            out2.flush();
            TransformerFactory transFactory = TransformerFactory.newInstance();
            Transformer transformer = transFactory.newTransformer();
            // Some options to do reformatting of the xml output; they don't seem to do anything to the mzid data.
            //transformer.setOutputProperty(OutputKeys.INDENT, "yes");
            //transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
            // Prevent automatic output of the xml declaration (also removes it if it is already present in the input)
            transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes");
            // Another way to make the xml declaration be on a separate line, but it also adds the string 'standalone="yes"' to the xml declaration
            //transformer.setOutputProperty(OutputKeys.STANDALONE, "yes");
            transformer.transform(new StreamSource(tempIn), new StreamResult(out2));
            out2.flush();
            out2.close();
            tempFile.deleteOnExit();
        } catch (IOException e) {
            System.out.println("Could not reprocess temp file \"" + tempFile + "\". Will not delete file (if file exists).");
            e.printStackTrace();
        } catch (TransformerException e) {
            e.printStackTrace();
        }
    }

    public MZIdentMLGen setEValueThreshold(float eValueThreshold) {
        this.eValueThreshold = eValueThreshold;
        return this;
    }

    public synchronized void addSpectrumIdentificationResults(List<MSGFPlusMatch> resultList) {
        for (MSGFPlusMatch mpMatch : resultList) {
            int specIndex = mpMatch.getSpecIndex();
            List<DatabaseMatch> matchList = mpMatch.getMatchList();
            if (matchList == null || matchList.size() == 0)
                continue;

            edu.ucsd.msjava.msutil.Spectrum spec = specAcc.getSpecMap().getSpectrumBySpecIndex(specIndex);
            String specID = spec.getID();
            float precursorMz = spec.getPrecursorPeak().getMz();

            SpectrumIdentificationResult sir = new SpectrumIdentificationResult();
            sir.setId(Constants.sirID + specIndex);
            sir.setSpectraData(spectraData);
            sir.setSpectrumID(specID);

            // add title
            String title = spec.getTitle();
            if (title != null) {
                CvParam cvParam = Constants.makeCvParam("MS:1000796", "spectrum title");
                cvParam.setValue(title);
                sir.getCvParam().add(cvParam);
            }

            // add scan number
            int scanNum = spec.getScanNum();
            if (scanNum >= 0) {
                CvParam cvParam = Constants.makeCvParam("MS:1001115", "scan number(s)");
                cvParam.setValue(String.valueOf(scanNum));
                sir.getCvParam().add(cvParam);
            }

            // add retention time
            float scanStartTime = spec.getRt();
            if (scanStartTime >= 0) {
                CvParam cvParam = Constants.makeCvParam("MS:1000016", "scan start time");
                cvParam.setValue(String.valueOf(scanStartTime));
                if (spec.getRtIsSeconds()) {
                    cvParam.setUnitCv(Constants.unitCV);
                    cvParam.setUnitAccession("UO:0000010");
                    cvParam.setUnitName("second");
                } else {
                    cvParam.setUnitCv(Constants.unitCV);
                    cvParam.setUnitAccession("UO:0000031");
                    cvParam.setUnitName("minute");
                }
                sir.getCvParam().add(cvParam);
            }


            int rank = 0;
            for (int i = matchList.size() - 1; i >= 0; --i) {
                ++rank;
                DatabaseMatch match = matchList.get(i);

                if (match.getDeNovoScore() < params.getMinDeNovoScore())
                    break;

//				int pepIndex = match.getIndex();	// Position of preAA
                int length = match.getLength();        // Peptide length + 2
                int charge = match.getCharge();

                float peptideMass = match.getPeptideMass();
//				float theoMass = peptideMass + (float)Composition.H2O;
//				float theoMz = theoMass/charge;
                float theoMz = (peptideMass + (float) Composition.H2O) / charge + (float) Composition.ChargeCarrierMass();

                int score = match.getScore();
                double specEValue = match.getSpecEValue();
                int numPeptides = sa.getNumDistinctPeptides(params.getEnzyme() == null ? length - 2 : length - 1);
                double eValue = specEValue * numPeptides;

                String specEValueStr;
                if (specEValue < Float.MIN_NORMAL)
                    specEValueStr = String.valueOf(specEValue);
                else
                    specEValueStr = String.valueOf((float) specEValue);

                String eValueStr;
                if (specEValue < Float.MIN_NORMAL)
                    eValueStr = String.valueOf(eValue);
                else
                    eValueStr = String.valueOf((float) eValue);

                SpectrumIdentificationItem sii = new SpectrumIdentificationItem();

                sii.setChargeState(charge);
                sii.setExperimentalMassToCharge(precursorMz);
                sii.setCalculatedMassToCharge((double) theoMz);

                Peptide pep = getPeptide(match);
                sii.setPeptide(pep);

                sii.setRank(rank);
                sii.setPassThreshold(eValue <= eValueThreshold);
                sii.setId(Constants.siiID + specIndex + "_" + rank);

                sii.getPeptideEvidenceRef().addAll(getPeptideEvidenceList(match, pep));

                List<CvParam> cvList = sii.getCvParam();
                List<UserParam> userList = sii.getUserParam();

//				ActivationMethod[] activationMethodArr = match.getActivationMethodArr();
//				if(activationMethodArr != null)
//				{
//					for(ActivationMethod actMethod : activationMethodArr)
//					{
//						CvParam fragMethodCV = actMethod.getCvParam();
//						cvList.add(fragMethodCV);
//					}
//				}

                CvParam rawScoreCV = Constants.makeCvParam("MS:1002049", "MS-GF:RawScore");
                rawScoreCV.setValue(String.valueOf(score));
                cvList.add(rawScoreCV);

                CvParam deNovoScoreCV = Constants.makeCvParam("MS:1002050", "MS-GF:DeNovoScore");
                deNovoScoreCV.setValue(String.valueOf(match.getDeNovoScore()));
                cvList.add(deNovoScoreCV);

                CvParam specEValueCV = Constants.makeCvParam("MS:1002052", "MS-GF:SpecEValue");
                specEValueCV.setValue(specEValueStr);
                cvList.add(specEValueCV);

                CvParam eValueCV = Constants.makeCvParam("MS:1002053", "MS-GF:EValue");
                eValueCV.setValue(eValueStr);
                cvList.add(eValueCV);

                if (match.getPSMQValue() != null) {
                    CvParam psmQValueCV = Constants.makeCvParam("MS:1002054", "MS-GF:QValue");
                    psmQValueCV.setValue(match.getPSMQValue().toString());
                    cvList.add(psmQValueCV);
                }

                if (match.getPepQValue() != null) {
                    CvParam pepQValueCV = Constants.makeCvParam("MS:1002055", "MS-GF:PepQValue");
                    pepQValueCV.setValue(match.getPepQValue().toString());
                    cvList.add(pepQValueCV);
                }

                UserParam isotopeErrorParam = Constants.makeUserParam("IsotopeError");
                float expMass = precursorMz * charge;
                float theoMass = theoMz * charge;
//				int isotopeError = NominalMass.toNominalMass(expMass) - NominalMass.toNominalMass(theoMass);
                int isotopeError = Math.round((expMass - theoMass) / (float) Composition.ISOTOPE);
                isotopeErrorParam.setValue(String.valueOf(isotopeError));
                userList.add(isotopeErrorParam);

                ActivationMethod[] activationMethodArr = match.getActivationMethodArr();
                if (activationMethodArr != null) {
                    StringBuffer actMethodStrBuf = new StringBuffer();
                    actMethodStrBuf.append(activationMethodArr[0]);
                    for (int j = 1; j < activationMethodArr.length; j++)
                        actMethodStrBuf.append("/" + activationMethodArr[j]);
                    UserParam dissociationMethodParam = Constants.makeUserParam("AssumedDissociationMethod");
                    dissociationMethodParam.setValue(actMethodStrBuf.toString());
                    userList.add(dissociationMethodParam);
                }

                if (match.getAdditionalFeatureList() != null) {
                    for (Pair<String, String> feature : match.getAdditionalFeatureList()) {
                        String name = feature.getFirst();
                        String value = feature.getSecond();
                        UserParam addParam = Constants.makeUserParam(name);
                        addParam.setValue(value);
                        userList.add(addParam);
                    }
                }

                sir.getSpectrumIdentificationItem().add(sii);
            }
            if (!sir.getSpectrumIdentificationItem().isEmpty())
                siList.getSpectrumIdentificationResult().add(sir);
        }
    }

    // index: peptide index in sarr
    // length: peptide length including pre/post
    // pepStr
    private Peptide getPeptide(DatabaseMatch match) {
//		String pvID = pre+pepStr+post;

        String pepStr = match.getPepSeq();

//		/////////
//		System.out.println("PepStr: " + pepStr);
//		///////////

        // pepStr is unique to a specific peptide/modification set
        Peptide mzidPeptide = pepMap.get(pepStr);

//		/////////////
//		edu.ucsd.msjava.msutil.Peptide testPep = aaSet.getPeptide(pepStr);
//		System.out.println(testPep);
//		System.out.println(pepStr+"->"+(mzidPeptide == null ? null : mzidPeptide.getId()));
        /////////////

        if (mzidPeptide == null) {
            // new peptide variant
            mzidPeptide = new Peptide();
            List<Modification> modList = mzidPeptide.getModification();
            edu.ucsd.msjava.msutil.Peptide peptide = aaSet.getPeptide(pepStr);
            StringBuffer unmodPepStr = new StringBuffer();
            StringBuffer modPepStr = new StringBuffer();
            int location = 1;
            for (edu.ucsd.msjava.msutil.AminoAcid aa : peptide) {
                Location loc;
                if (location == 1)
                    loc = Location.N_Term;
                else if (location == peptide.size())
                    loc = Location.C_Term;
                else
                    loc = Location.Anywhere;

                char residue = aa.getUnmodResidue();
                unmodPepStr.append(residue);

                List<Double> modMasses = new ArrayList<Double>();
                boolean hasNTermMod = false;
                List<Double> nTermMasses = new ArrayList<Double>();
                boolean hasCTermMod = false;
                List<Double> cTermMasses = new ArrayList<Double>();
                boolean modified = false;
                
                if (loc == Location.N_Term || loc == Location.C_Term) {
                    List<edu.ucsd.msjava.msutil.Modification> fixedTermMods = apcGen.getTerminalFixedModifications(aa.getUnmodResidue(), loc);
                    for (edu.ucsd.msjava.msutil.Modification fixedMod : fixedTermMods) {
                        Modification mod = new Modification();
                        if (loc == Location.N_Term)
                            mod.setLocation(location - 1);
                        else
                            mod.setLocation(location + 1);
                        mod.setMonoisotopicMassDelta(fixedMod.getAccurateMass());
                        mod.getCvParam().addAll(apcGen.getSearchModification(fixedMod).getCvParam());
                        modList.add(mod);
                    }
                }
                List<edu.ucsd.msjava.msutil.Modification> fixedMods = apcGen.getFixedModifications(aa.getUnmodResidue());
                if (fixedMods != null) {
                    for (edu.ucsd.msjava.msutil.Modification fixedMod : fixedMods) {
                        Modification mod = new Modification();
                        mod.setLocation(location);
                        mod.setMonoisotopicMassDelta(fixedMod.getAccurateMass());
                        mod.getCvParam().addAll(apcGen.getSearchModification(fixedMod).getCvParam());
                        modList.add(mod);
                    }
                }
                if (aa.isModified()) {
                    Modification mod = new Modification();
                    ModifiedAminoAcid modAA = (ModifiedAminoAcid) aa;
                    if (location == 1 && modAA.isNTermVariableMod()) {
                        mod.setLocation(location - 1);
                        nTermMasses.add(modAA.getModification().getAccurateMass());
                        hasNTermMod = true;
                    }
                    else if (location == peptide.size() && modAA.isCTermVariableMod()) {
                        mod.setLocation(location + 1);
                        cTermMasses.add(modAA.getModification().getAccurateMass());
                        hasCTermMod = true;
                    }
                    else {
                        mod.setLocation(location);
                        modMasses.add(modAA.getModification().getAccurateMass());
                        modified = true;
                    }
                    mod.setMonoisotopicMassDelta(modAA.getModification().getAccurateMass());

                    mod.getCvParam().addAll(apcGen.getSearchModification(modAA.getModification()).getCvParam());
                    modList.add(mod);

                    while (modAA.getTargetAA().isModified())    // aa has two (or more) modifications
                    {
                        Modification mod2 = new Modification();
                        modAA = (ModifiedAminoAcid) modAA.getTargetAA();
                        if (location == 1 && modAA.isNTermVariableMod()) {
                            mod2.setLocation(location - 1);
                            nTermMasses.add(modAA.getModification().getAccurateMass());
                            hasNTermMod = true;
                        }
                        else if (location == peptide.size() && modAA.isCTermVariableMod()) {
                            mod2.setLocation(location + 1);
                            cTermMasses.add(modAA.getModification().getAccurateMass());
                            hasCTermMod = true;
                        }
                        else {
                            mod2.setLocation(location);
                            modMasses.add(modAA.getModification().getAccurateMass());
                            modified = true;
                        }
                        mod2.setMonoisotopicMassDelta(modAA.getModification().getAccurateMass());

                        mod2.getCvParam().addAll(apcGen.getSearchModification(modAA.getModification()).getCvParam());
                        modList.add(mod2);
                    }
                }
                
                if (hasNTermMod) {
                    modPepStr.append("[");
                    Collections.sort(nTermMasses);
                    for (Double nTMass : nTermMasses) {
                        if (nTMass >= 0)
                            modPepStr.append("+");
                        modPepStr.append(Math.round(nTMass));
                    }
                }
                
                modPepStr.append(residue);
                if (modified) {
                    Collections.sort(modMasses);
                    for (Double mMass : modMasses) {
                        if (mMass >= 0)
                            modPepStr.append("+");
                        modPepStr.append(Math.round(mMass));
                    }
                }
                
                if (hasCTermMod) {
                    modPepStr.append("}");
                    Collections.sort(cTermMasses);
                    for (Double cTMass : cTermMasses) {
                        if (cTMass >= 0)
                            modPepStr.append("+");
                        modPepStr.append(Math.round(cTMass));
                    }
                }
                
                location++;
            }

            mzidPeptide.setPeptideSequence(unmodPepStr.toString());
            pepMap.put(pepStr, mzidPeptide);
            //mzidPeptide.setId(Constants.pepIDPrefix+pepMap.size());
            mzidPeptide.setId(Constants.pepIDPrefix+modPepStr.toString());
            peptideList.add(mzidPeptide);
        }

        return mzidPeptide;
    }

    public List<PeptideEvidenceRef> getPeptideEvidenceList(DatabaseMatch match, Peptide peptide) {
        SortedSet<Integer> indices = match.getIndices();
        int length = match.getLength();

        int startKey = indices.first();
//		if(match.isNTermMetCleaved())
//			++startKey;

        String pepIDNum = peptide.getId().substring(Constants.pepIDPrefix.length());
        String annotationKey = (match.isNTermMetCleaved() ? "M" : "") + startKey + "_" + pepIDNum;
        List<PeptideEvidenceRef> evRefList = evRefListMap.get(annotationKey);

        if (evRefList == null) {
            evRefList = new ArrayList<PeptideEvidenceRef>();

            CompactFastaSequence seq = sa.getSequence();
            for (int index : indices) {
                boolean isNTermMetCleaved;
                isNTermMetCleaved = match.isNTermMetCleaved() && sa.getSequence().getCharAt(index + 1) == 'M';

                PeptideEvidence pepEv = new PeptideEvidence();

                char pre = sa.getSequence().getCharAt(index);
                if (pre == '_') {
                    if (isNTermMetCleaved)
                        pre = 'M';
                    else
                        pre = '-';
                }
                char post;
                if (isNTermMetCleaved)
                    post = sa.getSequence().getCharAt(index + length);
                else
                    post = sa.getSequence().getCharAt(index + length - 1);
                if (post == '_')
                    post = '-';
                pepEv.setPre(String.valueOf(pre));
                pepEv.setPost(String.valueOf(post));

                int protStartIndex = (int) seq.getStartPosition(index);
                DBSequence dbSeq = getDBSequence(protStartIndex);
                pepEv.setDBSequence(dbSeq);
                pepEv.setPeptide(peptide);

                int start = index - protStartIndex + 1;
                if (isNTermMetCleaved)
                    ++start;

                int end = start + length - 2 - 1;
                pepEv.setStart(start);
                pepEv.setEnd(end);

//				String pepEvKey = "PepEv"+(index+1)+"_"+length;
                String pepEvKey = Constants.pepEvIDPrefix + (index + 1) + "_" + pepIDNum + "_" + start;
//				if(match.isNTermMetCleaved())
//					pepEvKey += "_N"+"_"+startKey;
                pepEv.setId(pepEvKey);

                pepEv.setIsDecoy(isDecoyMap.get(protStartIndex));

                if (pepEvMap.get(pepEvKey) != null) {
                    // Avoid duplicate peptide evidences
                    // currently only occurs when there are 2 results for a single peptide/mod/dbseq combo,
                    // one where match.isNTermMetCleaved() is true and sa.getSequence().getCharAt(index + 1) == 'M' is false,
                    // and another where match.isNTermMetCleaved() is false.
                    pepEv = pepEvMap.get(pepEvKey);
                }
                else {
                    pepEvMap.put(pepEvKey, pepEv);
                }

                peptideEvidenceList.add(pepEv);

                PeptideEvidenceRef pepEvRef = new PeptideEvidenceRef();
                pepEvRef.setPeptideEvidence(pepEv);
                evRefList.add(pepEvRef);
            }
            evRefListMap.put(annotationKey, evRefList);
        }

        return evRefList;
    }

    public DBSequence getDBSequence(int protStartIndex) {
        DBSequence dbSeq = dbSeqMap.get(protStartIndex);
        if (dbSeq == null) {
            dbSeq = new DBSequence();

            CompactFastaSequence seq = sa.getSequence();
            String annotation = seq.getAnnotation(protStartIndex);
            String proteinSeq = seq.getMatchingEntry(protStartIndex);
            String accession = annotation.split("\\s+")[0];

            dbSeq.setLength(proteinSeq.length());
            dbSeq.setSearchDatabase(searchDatabase);
            dbSeq.setAccession(accession);
            dbSeq.setId("DBSeq" + (protStartIndex + 1));

            boolean isDecoy = accession.startsWith(MSGFPlus.DECOY_PROTEIN_PREFIX);
            if (!isDecoy) {
                CvParam protDescCV = Constants.makeCvParam("MS:1001088", "protein description");
                protDescCV.setValue(annotation);
                dbSeq.getCvParam().add(protDescCV);
            }

            this.dbSequenceList.add(dbSeq);
            dbSeqMap.put(protStartIndex, dbSeq);
            isDecoyMap.put(protStartIndex, isDecoy);
        }

        return dbSeq;
    }

    private void generateCvList() {
        cvList = new CvList();
        List<Cv> localCvList = cvList.getCv();
        localCvList.add(Constants.psiCV);
        localCvList.add(Constants.unimodCV);
        localCvList.add(Constants.unitCV);
    }

    private void generateAnalysisSoftwareList() {
        analysisSoftwareList = new AnalysisSoftwareList();
        List<AnalysisSoftware> analysisSoftwares = analysisSoftwareList.getAnalysisSoftware();
        analysisSoftwares.add(Constants.msgfPlus);
    }

    private void initSequenceCollection() {
        sequenceCollection = new SequenceCollection();
        dbSequenceList = sequenceCollection.getDBSequence();
        peptideList = sequenceCollection.getPeptide();
        peptideEvidenceList = sequenceCollection.getPeptideEvidence();
    }

    private void generateAnalysisCollection() {

        analysisCollection = new AnalysisCollection();

        List<SpectrumIdentification> specIdentList = analysisCollection.getSpectrumIdentification();

        SpectrumIdentification specIdent = new SpectrumIdentification();
        specIdent.setId(Constants.specIdentID);
        specIdent.setSpectrumIdentificationList(siList);
        specIdent.setSpectrumIdentificationProtocol(siProtocol);

        List<InputSpectra> inputSpecList = specIdent.getInputSpectra();
        InputSpectra inputSpec = new InputSpectra();
        inputSpec.setSpectraData(spectraData);
        inputSpecList.add(inputSpec);


        List<SearchDatabaseRef> searchDBRefList = specIdent.getSearchDatabaseRef();
        SearchDatabaseRef searchDBRef = new SearchDatabaseRef();
        searchDBRef.setSearchDatabase(searchDatabase);
        searchDBRefList.add(searchDBRef);

        specIdentList.add(specIdent);
    }

    private void generateAnalysisProtocolCollection() {
        apcGen = new AnalysisProtocolCollectionGen(params, aaSet);
        analysisProtocolCollection = apcGen.getAnalysisProtocolCollection();
        siProtocol = apcGen.getSpectrumIdentificationProtocol();
    }

    private void initDataCollection() {
        dataCollection = new DataCollection();

        // Inputs
        Inputs inputs = new Inputs();
        // source file: skip

        // search database
        searchDatabase = new SearchDatabase();
        searchDatabase.setId(Constants.searchDBID);
        searchDatabase.setNumDatabaseSequences((long) sa.getSequence().getNumProteins());
        searchDatabase.setLocation(params.getDatabaseFile().getAbsolutePath());

        UserParam param = new UserParam();
        param.setName(params.getDatabaseFile().getName());
        Param tempParam = new Param();
        tempParam.setParam(param);
        searchDatabase.setDatabaseName(tempParam);

        FileFormat ffDB = new FileFormat();
        ffDB.setCvParam(Constants.makeCvParam("MS:1001348", "FASTA format"));
        searchDatabase.setFileFormat(ffDB);

        // for decoy
        if (params.useTDA()) {
            searchDatabase.getCvParam().add(Constants.makeCvParam("MS:1001197", "DB composition target+decoy"));
            CvParam decoyAccCV = Constants.makeCvParam("MS:1001283", "decoy DB accession regexp");
            decoyAccCV.setValue("^" + MSGFPlus.DECOY_PROTEIN_PREFIX);
            searchDatabase.getCvParam().add(decoyAccCV);
            searchDatabase.getCvParam().add(Constants.makeCvParam("MS:1001195", "decoy DB type reverse"));
        }

        inputs.getSearchDatabase().add(searchDatabase);

        // spectra data
        spectraData = new SpectraData();
        spectraData.setId(Constants.spectraDataID);

        File specFile = params.getDBSearchIOList().get(ioIndex).getSpecFile();
        spectraData.setLocation(specFile.getAbsolutePath());
        spectraData.setName(specFile.getName());

        // spectrum file format, TODO: add _dta.txt to the PSI CV
        SpecFileFormat specFileFormat = params.getDBSearchIOList().get(ioIndex).getSpecFileFormat();
        FileFormat ffSpec = new FileFormat();
        ffSpec.setCvParam(Constants.makeCvParam(specFileFormat.getPSIAccession(), specFileFormat.getPSIName()));
        spectraData.setFileFormat(ffSpec);

        SpectrumIDFormat sidFormat = new SpectrumIDFormat();
        sidFormat.setCvParam(specAcc.getSpectrumIDFormatCvParam());
        spectraData.setSpectrumIDFormat(sidFormat);

        inputs.getSpectraData().add(spectraData);
        dataCollection.setInputs(inputs);

        // AnalysisData
        AnalysisData analysisData = new AnalysisData();
        dataCollection.setAnalysisData(analysisData);

        siList = new SpectrumIdentificationList();
        siList.setId(Constants.siListID);

        FragmentationTable fragTable = new FragmentationTable();
        List<Measure> measureList = fragTable.getMeasure();
        Measure mzMeasure = new Measure();
        mzMeasure.setId(Constants.measureMzID);
        List<CvParam> cvParamList = mzMeasure.getCvParam();
        cvParamList.add(Constants.makeCvParam("MS:1001225", "product ion m/z", Constants.psiCV, "MS:1000040", "m/z", Constants.psiCV));
        measureList.add(mzMeasure);

//		Measure intMeasure = new Measure();
//		intMeasure.setId(Constants.measureIntID);
//		cvParamList = intMeasure.getCvParam();
//		cvParamList.add(Constants.makeCvParam("MS:1001226","product ion intensity", Constants.psiCV,"MS:1000131", "number of counts", Constants.psiCV));
//		measureList.add(intMeasure);

//		Measure errorMeasure = new Measure();
//		errorMeasure.setId(Constants.measureErrorID);
//		cvParamList = errorMeasure.getCvParam();
//		cvParamList.add(Constants.makeCvParam("MS:1001227","product ion m/z error", Constants.psiCV,"MS:1000040","m/z", Constants.psiCV));
//		measureList.add(errorMeasure);

        siList.setFragmentationTable(fragTable);
        analysisData.getSpectrumIdentificationList().add(siList);
    }


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
