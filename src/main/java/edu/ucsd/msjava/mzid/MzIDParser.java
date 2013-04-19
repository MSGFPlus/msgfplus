package edu.ucsd.msjava.mzid;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.mzml.MzMLAdapter;

import uk.ac.ebi.jmzidml.model.mzidml.AnalysisData;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.DataCollection;
import uk.ac.ebi.jmzidml.model.mzidml.FileFormat;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidenceRef;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectraData;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationProtocol;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.model.mzidml.Tolerance;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

public class MzIDParser {
	private final MzIdentMLUnmarshaller unmarshaller;
	private final boolean showDecoy;
	private boolean doNotShowQValue;
	private final boolean unrollResults;
	private final boolean showMolecularFormula;
	
	private boolean isPrecursorTolerancePPM;
	private Map<String, Peptide> pepMap;		// Peptide ref -> Peptide
	private Map<String, DBSequence> dbSeqMap;		// DBSequenhce ref -> DBSequence
	private Map<String, PeptideEvidence> pepEvMap;	// PeptideEvidence ref -> PeptideEvidence

	public MzIDParser(File mzIDFile)
	{
		this(mzIDFile, false, false, false, false);
	}

	public MzIDParser(File mzIDFile, boolean showDecoy, boolean doNotShowQValue, boolean unrollResults, boolean showMolecularFormula)
	{
		unmarshaller = new MzIdentMLUnmarshaller(mzIDFile);
		this.showDecoy = showDecoy;
		this.doNotShowQValue = doNotShowQValue;
		this.unrollResults = unrollResults;
		this.showMolecularFormula = showMolecularFormula;
	}
	
	public void writeToTSVFile(File outputFile)
	{
		PrintStream out = null;
		if(outputFile != null)
			try {
				out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		else
			out = System.out;
		
		writeToTSVFile(out);
		
		if(out != System.out)
			out.close();
	}
	
	public void writeToTSVFile(PrintStream out)
	{
		unmarshallSequenceCollection();
		unmarshallAnalysisProtocolCollection();
		
        DataCollection dc =  unmarshaller.unmarshal(DataCollection.class);
        
        // get spectrum file
        boolean isMgf = false;
        Map<String, String> specFileNameMap = new HashMap<String, String>();
        Inputs inputs = dc.getInputs();
        for(SpectraData sd : inputs.getSpectraData())
        {
        	String specFileName = new File(sd.getLocation()).getName();
        	specFileNameMap.put(sd.getId(), specFileName);
        	FileFormat ff = sd.getFileFormat();
        	if(ff.getCvParam().getAccession().equals("MS:1001062"))
        		isMgf = true;
        }

		String header = 
				"#SpecFile" +
				"\tSpecID" +
				"\tScanNum" +
				(isMgf ? "\tTitle" : "") +
				"\tFragMethod"
				+"\tPrecursor"
				+"\tIsotopeError"
				+"\tPrecursorError("
				+ (isPrecursorTolerancePPM ? "ppm" : "Da")
				+")" +
				"\tCharge" +
				"\tPeptide" +
				(showMolecularFormula ? "\tFormula" : "") +
				"\tProtein" +
				"\tDeNovoScore" +
				"\tMSGFScore" +
				"\tSpecEValue" +
				"\tEValue" +
				(!this.doNotShowQValue ? "\tQValue\tPepQValue" : "");
		out.println(header);
        
        AnalysisData ad = dc.getAnalysisData();

        // Get the list of SpectrumIdentification elements
        List<SpectrumIdentificationList> sil = ad.getSpectrumIdentificationList();

        for (SpectrumIdentificationList sIdentList : sil) {
             for (SpectrumIdentificationResult sir
                     : sIdentList.getSpectrumIdentificationResult()) {

                 String specFileName = specFileNameMap.get(sir.getSpectraDataRef());
            	 String specID = sir.getSpectrumID();
            	 Map<String, CvParam> sirCvParamMap = getCvParamMap(sir.getCvParam());
            	 
            	 String scanNum = "-1";
            	 CvParam scanNumParam = sirCvParamMap.get("MS:1001115");
            	 if(scanNumParam != null)
            		 scanNum = scanNumParam.getValue();

            	 String title = "N/A";
            	 CvParam titleParam = sirCvParamMap.get("MS:1000796");
            	 if(titleParam != null)
            		 title = titleParam.getValue();
            	 
                 for (SpectrumIdentificationItem sii
                      : sir.getSpectrumIdentificationItem()) {

                     Double calculatedMassToCharge =  sii.getCalculatedMassToCharge();
                     Double experimentalMassToCharge = sii.getExperimentalMassToCharge();
                     int charge = sii.getChargeState();

                     Map<String, CvParam> cvParamMap = getCvParamMap(sii.getCvParam());
                     CvParam cvParam;
                     
//                     String fragMethod = null;
//                     for(ActivationMethod actMethod : ActivationMethod.getAllRegisteredActivationMethods())
//                     {
//                    	 if((cvParam = cvParamMap.get(actMethod.getPSICVAccession())) != null)
//                    	 {
//                    		 fragMethod = actMethod.getName();
//                    		 break;
//                    	 }
//                     }

                     
                     String deNovoScore =  (cvParam = cvParamMap.get("MS:1002050")) == null ? "" : cvParam.getValue();
                     String rawScore =  (cvParam = cvParamMap.get("MS:1002049")) == null ? "" : cvParam.getValue();
                     String specEValue =  (cvParam = cvParamMap.get("MS:1002052")) == null ? "" : cvParam.getValue();
                     String eValue =  (cvParam = cvParamMap.get("MS:1002053")) == null ? "" : cvParam.getValue();
                     String psmQValue =  (cvParam = cvParamMap.get("MS:1002054")) == null ? "" : cvParam.getValue();
                     String pepQValue =  (cvParam = cvParamMap.get("MS:1002055")) == null ? "" : cvParam.getValue();

                     Map<String, UserParam> userParamMap = getUserParamMap(sii.getUserParam());
                     UserParam userParam;
                     String fragMethod =  (userParam = userParamMap.get("AssumedDissociationMethod")) == null ? null : userParam.getValue();
                     Integer isotopeError =  (userParam = userParamMap.get("IsotopeError")) == null ? null : Integer.parseInt(userParam.getValue());
                     
                     double adjustedExpMz = experimentalMassToCharge-Composition.ISOTOPE*isotopeError/charge;
                     double precursorError = adjustedExpMz-calculatedMassToCharge;
                     if(isPrecursorTolerancePPM)
                    	 precursorError = precursorError/calculatedMassToCharge*1e6;
                     
                     Peptide peptide = pepMap.get(sii.getPeptideRef());
                     String peptideSeq = getPeptideSeq(peptide);
                     String molecularFormula = null;
                     if(showMolecularFormula)
                    	 molecularFormula = getMolecularFormula(peptide);

                	 HashSet<String> proteinSet = new HashSet<String>();
                     if(this.unrollResults)
                     {
                         for(PeptideEvidenceRef pepEvRef : sii.getPeptideEvidenceRef())
                         {
                        	 PeptideEvidence pepEv = pepEvMap.get(pepEvRef.getPeptideEvidenceRef());
                        	 
                        	 boolean isDecoy = pepEv.isIsDecoy();
                        	 if(isDecoy && !this.showDecoy)
                        		 continue;
                        	 
                        	 String pre = pepEv.getPre();
                        	 String post = pepEv.getPost();
                        	 
                        	 DBSequence dbSeq = dbSeqMap.get(pepEv.getDBSequenceRef());
                        	 String protein = dbSeq.getAccession();
                        	 if(proteinSet.add(pre+protein+post))
                        	 {
                                 out.print(specFileName
                                		 +"\t"+specID
                                		 +"\t"+scanNum
                                		 +(!isMgf ? "" : "\t"+title)
                                		 +"\t"+fragMethod
                                		 +"\t"+experimentalMassToCharge.floatValue()
                                		 +"\t"+isotopeError
                                		 +"\t"+(float)precursorError
                                		 +"\t"+charge
                                		 +"\t"+pre+"."+peptideSeq+"."+post
                                		 +(molecularFormula == null ? "" : "\t"+molecularFormula)
                                		 +"\t"+protein
                                		 +"\t"+deNovoScore
                                		 +"\t"+rawScore
                                		 +"\t"+specEValue
                                		 +"\t"+eValue
                                		 );
                                 if(!this.doNotShowQValue)
                                	 out.print("\t"+psmQValue+"\t"+pepQValue);
                                 out.println();
                        	 }
                         }                  	 
                     }
                     else
                     {
                    	 StringBuffer proteinBuf = new StringBuffer();
                    	 boolean isAllDecoy = true;
                         for(PeptideEvidenceRef pepEvRef : sii.getPeptideEvidenceRef())
                         {
                        	 PeptideEvidence pepEv = pepEvMap.get(pepEvRef.getPeptideEvidenceRef());
                        	 
                        	 boolean isDecoy = pepEv.isIsDecoy();
                        	 if(isDecoy && !this.showDecoy)
                        	 {
                        		 continue;
                        	 }
                        	 
                        	 isAllDecoy = false;
                        	 String pre = pepEv.getPre();
                        	 String post = pepEv.getPost();
                        	 
                        	 DBSequence dbSeq = dbSeqMap.get(pepEv.getDBSequenceRef());
                        	 String protein = dbSeq.getAccession();
                        	 
                        	 if(proteinSet.add(pre+protein+post))
                        	 {
                        		 if(proteinBuf.length() != 0)
                        			 proteinBuf.append(";");
                        		 proteinBuf.append(protein+"(pre="+pre+",post="+post+")");
                        	 }
                         }
                         
                         if(!isAllDecoy)
                         {
                             out.print(specFileName
                            		 +"\t"+specID
                            		 +"\t"+scanNum
                            		 +(!isMgf ? "" : "\t"+title)
                            		 +"\t"+fragMethod
                            		 +"\t"+experimentalMassToCharge.floatValue()
                            		 +"\t"+isotopeError
                            		 +"\t"+(float)precursorError
                            		 +"\t"+charge
                            		 +"\t"+peptideSeq
                            		 +(molecularFormula == null ? "" : "\t"+molecularFormula)
                            		 +"\t"+proteinBuf.toString()
                            		 +"\t"+deNovoScore
                            		 +"\t"+rawScore
                            		 +"\t"+specEValue
                            		 +"\t"+eValue
                            		 );
	                         if(!this.doNotShowQValue)
	                        	 out.print("\t"+psmQValue+"\t"+pepQValue);
	                         out.println();
                         }
                     }

                 } // end spectrum identification item
             } // end spectrum identification results
        }
	}
	
	private Map<String,CvParam> getCvParamMap(List<CvParam> paramList)
	{
		Map<String,CvParam> paramMap = new HashMap<String, CvParam>();
		
        for(CvParam param : paramList)
        	paramMap.put(param.getAccession(), param);
        
		return paramMap;
	}

	private Map<String,UserParam> getUserParamMap(List<UserParam> paramList)
	{
		Map<String,UserParam> paramMap = new HashMap<String, UserParam>();
		
        for(UserParam param : paramList)
        	paramMap.put(param.getName(), param);
        
		return paramMap;
	}

	private void unmarshallSequenceCollection()
	{
        SequenceCollection sc = unmarshaller.unmarshal(SequenceCollection.class);
        
        dbSeqMap = new HashMap<String, DBSequence>();
        for(DBSequence dbSeq : sc.getDBSequence())
        	dbSeqMap.put(dbSeq.getId(), dbSeq);
        
        pepMap = new HashMap<String, Peptide>();
        for(Peptide peptide : sc.getPeptide())
        	pepMap.put(peptide.getId(), peptide);
        
        pepEvMap = new HashMap<String, PeptideEvidence>();
        for(PeptideEvidence pepEv : sc.getPeptideEvidence())
        	pepEvMap.put(pepEv.getId(), pepEv);
	}
	
	private void unmarshallAnalysisProtocolCollection()
	{
		AnalysisProtocolCollection apc = unmarshaller.unmarshal(AnalysisProtocolCollection.class);
		SpectrumIdentificationProtocol sip = apc.getSpectrumIdentificationProtocol().get(0);
		
		isPrecursorTolerancePPM = false;
		Tolerance parentTolerance = sip.getParentTolerance();
		for(CvParam param : parentTolerance.getCvParam())
		{
			if(param.getAccession().equals("MS:1001412"))
			{
				if(param.getUnitName().equals("parts per million"))
				{
					isPrecursorTolerancePPM = true;
					break;
				}
			}
		}

		if(!doNotShowQValue)
		{
			for(UserParam param : sip.getAdditionalSearchParams().getUserParam())
			{
				if(param.getName().equals("TargetDecoyApproach"))
				{
					if(param.getValue().equals("false"))
					{
						doNotShowQValue = true;
					}
					break;
				}
			}
		}
	}
	
	private static String getPeptideSeq(Peptide peptide)
	{
		String unmodPepSeq = peptide.getPeptideSequence();
		
		String[] modArr = new String[unmodPepSeq.length()+2];
        // Modification
        for(Modification mod : peptide.getModification())
        {
        	double modMass = mod.getMonoisotopicMassDelta();
        	String massStr;
    		if(modMass >= 0)
    			massStr = "+" + String.format("%.3f", modMass);
    		else
    			massStr = String.format("%.3f", modMass);
    		if(modArr[mod.getLocation()] == null)
    			modArr[mod.getLocation()] = massStr;
    		else
    			modArr[mod.getLocation()] += massStr;
        }

		StringBuffer buf = new StringBuffer();
		if(modArr[0] != null)
			buf.append(modArr[0]);
        for(int i=0; i<unmodPepSeq.length(); i++)
        {
        	buf.append(unmodPepSeq.charAt(i));
        	if(modArr[i+1] != null)
        		buf.append(modArr[i+1]);
        }
        if(modArr[modArr.length-1] != null)
        	buf.append(modArr[modArr.length-1]);
		
		return buf.toString();
	}
	
	private static String getMolecularFormula(Peptide peptide)
	{
		AminoAcidSet stdAASet = AminoAcidSet.getStandardAminoAcidSet();
		
		String unmodPepSeq = peptide.getPeptideSequence();
		UnimodComposition composition = new UnimodComposition();
		for(int i=0; i<unmodPepSeq.length(); i++)
		{
			char residue = unmodPepSeq.charAt(i);
			composition.add(stdAASet.getAminoAcid(residue).getComposition());
		}
		
        // Modification
        for(Modification mod : peptide.getModification())
        {
        	boolean hasComposition = false;
        	for(CvParam cvParam : mod.getCvParam())
        	{
        		String accession = cvParam.getAccession();
    			String deltaComposition = Unimod.getUnimod().getDeltaComposition(accession);
    			if(deltaComposition != null)	// correct unimod accession number
    			{
    				composition.add(deltaComposition);
    				hasComposition = true;
    				break;
    			}
        	}
        	if(!hasComposition)
        		composition.add(mod.getMonoisotopicMassDelta());
        }
        
        composition.add("H", 2);	// add H2O
        composition.add("O", 1);
        
        return composition.toString();
	}


	public static void main(String argv[]) throws Exception
	{
		long time = System.currentTimeMillis();
		MzMLAdapter.turnOffLogs();
//		File mzidFile = new File(System.getProperty("user.home")+"/Research/Data/QCShew/test.mzid");
//		File outputFile = null;
		
		File mzidFile = new File(System.getProperty("user.home")+"/Research/Data/QCShew/test.mzid");
		File outputFile = new File(System.getProperty("user.home")+"/Research/Data/QCShew/test.tsv");
		MzIDParser parser = new MzIDParser(mzidFile);
		parser.writeToTSVFile(outputFile);
		
		System.out.println("Elapsed time: " + (System.currentTimeMillis()-time)/1000f);
	}	

}
