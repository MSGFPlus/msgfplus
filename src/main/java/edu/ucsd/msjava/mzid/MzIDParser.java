package edu.ucsd.msjava.mzid;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.mzml.MzMLAdapter;

import uk.ac.ebi.jmzidml.model.mzidml.AnalysisData;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.DataCollection;
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
	
	private boolean isPrecursorTolerancePPM;
	private boolean isTDA;
	private Map<String, Peptide> pepMap;		// Peptide ref -> Peptide
	private Map<String, DBSequence> dbSeqMap;		// DBSequenhce ref -> DBSequence
	private Map<String, PeptideEvidence> pepEvMap;	// PeptideEvidence ref -> PeptideEvidence

	public MzIDParser(File mzIDFile)
	{
		unmarshaller = new MzIdentMLUnmarshaller(mzIDFile);
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
		
		unmarshallSequenceCollection();
		unmarshallAnalysisProtocolCollection();
		
		String header = 
				"#SpecFile" +
				"\tSpecID" +
				"\tScanNum" +
				"\tFragMethod"
				+"\tPrecursor" +
				"\tPrecursorError("
				+ (isPrecursorTolerancePPM ? "ppm" : "Da")
				+")" +
				"\tCharge" +
				"\tPeptide" +
				"\tProtein" +
				"\tDeNovoScore" +
				"\tMSGFScore" +
				"\tSpecEValue" +
				"\tEValue" +
				(this.isTDA ? "\tQValue\tPepQValue" : "");
		out.println(header);
		
        DataCollection dc =  unmarshaller.unmarshal(DataCollection.class);
        
        // get spectrum file
        Map<String, String> specFileNameMap = new HashMap<String, String>();
        Inputs inputs = dc.getInputs();
        for(SpectraData sd : inputs.getSpectraData())
        {
        	String specFileName = new File(sd.getLocation()).getName();
        	specFileNameMap.put(sd.getId(), specFileName);
        }
        
        AnalysisData ad = dc.getAnalysisData();

        // Get the list of SpectrumIdentification elements
        List<SpectrumIdentificationList> sil = ad.getSpectrumIdentificationList();

        for (SpectrumIdentificationList sIdentList : sil) {
             for (SpectrumIdentificationResult sir
                     : sIdentList.getSpectrumIdentificationResult()) {

                 String specFileName = specFileNameMap.get(sir.getSpectraDataRef());
            	 String specID = sir.getSpectrumID();
            	 Map<String, CvParam> sirCvParamMap = getCvParamMap(sir.getCvParam());
            	 
            	 String scanNum = "";
            	 CvParam scanNumParam = sirCvParamMap.get("MS:1001115");
            	 if(scanNumParam != null)
            		 scanNum = scanNumParam.getValue();
            	 
                 for (SpectrumIdentificationItem sii
                      : sir.getSpectrumIdentificationItem()) {

                     Double calculatedMassToCharge =  sii.getCalculatedMassToCharge() + Composition.H;
                     Double experimentalMassToCharge = sii.getExperimentalMassToCharge() + Composition.H;
                     int charge = sii.getChargeState();

                     Map<String, CvParam> cvParamMap = getCvParamMap(sii.getCvParam());
                     CvParam cvParam;
                     
                     String fragMethod = null;
                     for(ActivationMethod actMethod : ActivationMethod.getAllRegisteredActivationMethods())
                     {
                    	 if((cvParam = cvParamMap.get(actMethod.getPSICVAccession())) != null)
                    	 {
                    		 fragMethod = actMethod.getName();
                    		 break;
                    	 }
                     }
                     
                     String deNovoScore =  (cvParam = cvParamMap.get("MS:1002050")) == null ? "" : cvParam.getValue();
                     String rawScore =  (cvParam = cvParamMap.get("MS:1002049")) == null ? "" : cvParam.getValue();
                     String specEValue =  (cvParam = cvParamMap.get("MS:1002052")) == null ? "" : cvParam.getValue();
                     String eValue =  (cvParam = cvParamMap.get("MS:1002053")) == null ? "" : cvParam.getValue();
                     String psmQValue =  (cvParam = cvParamMap.get("MS:1002054")) == null ? "" : cvParam.getValue();
                     String pepQValue =  (cvParam = cvParamMap.get("MS:1002055")) == null ? "" : cvParam.getValue();

                     Map<String, UserParam> userParamMap = getUserParamMap(sii.getUserParam());
                     UserParam userParam;
                     Integer isotopeError =  (userParam = userParamMap.get("IsotopeError")) == null ? null : Integer.parseInt(userParam.getValue());
                     
                     double theoreticalMz = experimentalMassToCharge-Composition.ISOTOPE*isotopeError/charge;
                     double precursorError = calculatedMassToCharge-theoreticalMz;
                     if(isPrecursorTolerancePPM)
                    	 precursorError = precursorError/theoreticalMz*1e6;
                     
                     Peptide peptide = pepMap.get(sii.getPeptideRef());
                     String peptideSeq = getPeptideSeq(peptide);

                     for(PeptideEvidenceRef pepEvRef : sii.getPeptideEvidenceRef())
                     {
                    	 PeptideEvidence pepEv = pepEvMap.get(pepEvRef.getPeptideEvidenceRef());
                    	 String pre = pepEv.getPre();
                    	 String post = pepEv.getPost();
                    	 
                    	 DBSequence dbSeq = dbSeqMap.get(pepEv.getDBSequenceRef());
                    	 String protein = dbSeq.getAccession();

                         out.print(specFileName
                        		 +"\t"+specID
                        		 +"\t"+scanNum
                        		 +"\t"+fragMethod
                        		 +"\t"+calculatedMassToCharge.floatValue()
                        		 +"\t"+(float)precursorError
                        		 +"\t"+charge
                        		 +"\t"+pre+"."+peptideSeq+"."+post
                        		 +"\t"+protein
                        		 +"\t"+deNovoScore
                        		 +"\t"+rawScore
                        		 +"\t"+specEValue
                        		 +"\t"+eValue
                        		 );
                         if(this.isTDA)
                        	 out.print("\t"+psmQValue+"\t"+pepQValue);
                         out.println();
                    	 
                     }
                 } // end spectrum identification item
             } // end spectrum identification results
        }
        
        if(out != System.out)
        	out.close();
	}
	
	public Map<String,CvParam> getCvParamMap(List<CvParam> paramList)
	{
		Map<String,CvParam> paramMap = new HashMap<String, CvParam>();
		
        for(CvParam param : paramList)
        	paramMap.put(param.getAccession(), param);
        
		return paramMap;
	}

	public Map<String,UserParam> getUserParamMap(List<UserParam> paramList)
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

		isTDA = false;
		for(UserParam param : sip.getAdditionalSearchParams().getUserParam())
		{
			if(param.getName().equals("TargetDecoyApproach"))
			{
				if(param.getValue().equals("true"))
				{
					isTDA = true;
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
	


	public static void main(String argv[]) throws Exception
	{
		long time = System.currentTimeMillis();
		MzMLAdapter.turnOffLogs();
//		File mzidFile = new File(System.getProperty("user.home")+"/Research/Data/QCShew/test.mzid");
//		File outputFile = null;
		
		File mzidFile = new File(System.getProperty("user.home")+"/Research/Data/QCShew/QC_Shew_MSGFPlus_N3.mzid");
		File outputFile = new File(System.getProperty("user.home")+"/Research/Data/QCShew/QC_Shew_MSGFPlus_N3.tsv");
		MzIDParser parser = new MzIDParser(mzidFile);
		parser.writeToTSVFile(outputFile);
		
		System.out.println("Elapsed time: " + (System.currentTimeMillis()-time)/1000f);
	}	

}
