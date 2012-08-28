package edu.ucsd.msjava.mzid;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.mzml.MzMLAdapter;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisData;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.DataCollection;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

public class MzIDParser {
	private final File mzIDFile;
	private final MzIdentMLUnmarshaller unmarshaller;
	
	public MzIDParser(File mzIDFile)
	{
		this.mzIDFile = mzIDFile;
		unmarshaller = new MzIdentMLUnmarshaller(mzIDFile);
	}
	
	public void test()
	{
        DataCollection dc =  unmarshaller.unmarshal(DataCollection.class);
        AnalysisData ad = dc.getAnalysisData();

        // Get the list of SpectrumIdentification elements
        List<SpectrumIdentificationList> sil = ad.getSpectrumIdentificationList();

        for (SpectrumIdentificationList sIdentList : sil) {
             for (SpectrumIdentificationResult spectrumIdentResult
                     : sIdentList.getSpectrumIdentificationResult()) {

                 // Get the name of SpectrumIdentificationResult
                 String spectrumID =  spectrumIdentResult.getSpectrumID();

                 for (SpectrumIdentificationItem spectrumIdentItem
                      : spectrumIdentResult.getSpectrumIdentificationItem()) {

                     // Get the following information for SpectrumIdentificationItem element
                     String spectrumIdItem = spectrumIdentItem.getId();
                     Double calculatedMassToCharge =  spectrumIdentItem.getCalculatedMassToCharge();
                     Double experimentalMassToCharge = spectrumIdentItem.getExperimentalMassToCharge();
                     int rank = spectrumIdentItem.getRank();
                     int charge = spectrumIdentItem.getChargeState();

                     System.out.println("Spectrum Identification ID = " + spectrumIdItem);
                     System.out.println("Calculated Mass/Charge = " + calculatedMassToCharge);
                     System.out.println("Experimental Mass/Charge = " + experimentalMassToCharge);
                     System.out.println("Search rank = " + rank);
                     System.out.println("Charge = " + charge);

                     // If the auto-resolve mechanism is activated for SpectrumIdentificationItem
                     // then automatically resolve the Peptide Object
                     if (MzIdentMLElement.SpectrumIdentificationItem.isAutoRefResolving()
                             && spectrumIdentItem.getPeptideRef() != null) {
                          Peptide peptide = spectrumIdentItem.getPeptide();
                          String peptideSequence = peptide.getPeptideSequence();

                         System.out.println("Pepetide Sequence = " + peptideSequence);
                     }

                     System.out.println("\n");

                 } // end spectrum identification item
             } // end spectrum identification results
        }		
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

	Map<String, Peptide> pepMap;		// Peptide ref -> Peptide
	Map<String, DBSequence> dbSeqMap;		// DBSequenhce ref -> DBSequence
	Map<String, PeptideEvidence> pepEvMap;	// PeptideEvidence ref -> PeptideEvidence
	
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
	
	public static String getPeptideSeq(Peptide peptide)
	{
		StringBuffer buf = new StringBuffer();
		
		
		return buf.toString();
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
		
		String header = 
				"#SpecFile" +
				"\tSpecID" +
				"\tFragMethod"
				+"\tPrecursor" +
				"\tPrecursorError("
//				+(rightParentMassTolerance.isTolerancePPM() ? "ppm" : "Da")
				+"Da"
				+")" +
				"\tCharge" +
				"\tPeptide" +
				"\tProtein" +
				"\tDeNovoScore" +
				"\tMSGFScore" +
				"\tSpecEValue" +
				"\tEValue";
		out.println(header);
		
		unmarshallSequenceCollection();
        DataCollection dc =  unmarshaller.unmarshal(DataCollection.class);
        
        AnalysisData ad = dc.getAnalysisData();

        // Get the list of SpectrumIdentification elements
        List<SpectrumIdentificationList> sil = ad.getSpectrumIdentificationList();

        String specFile = "null";
        
        for (SpectrumIdentificationList sIdentList : sil) {
             for (SpectrumIdentificationResult sir
                     : sIdentList.getSpectrumIdentificationResult()) {

            	 String specID = sir.getSpectrumID();
                 for (SpectrumIdentificationItem sii
                      : sir.getSpectrumIdentificationItem()) {

                     Double calculatedMassToCharge =  sii.getCalculatedMassToCharge() + Composition.H;
                     Double experimentalMassToCharge = sii.getExperimentalMassToCharge() + Composition.H;
                     int charge = sii.getChargeState();

                     Peptide peptide = pepMap.get(sii.getPeptideRef());
                     String peptideSeq = peptide.getPeptideSequence();
                     
                     sii.getPeptideEvidenceRef();
                     
                     // Modification
                     peptide.getModification();
                     
                     Map<String, CvParam> cvParamMap = getCvParamMap(sii.getCvParam());
                     CvParam cvParam;
                     String deNovoScore =  (cvParam = cvParamMap.get("MS:1002050")) == null ? "null" : cvParam.getValue();
                     String rawScore =  (cvParam = cvParamMap.get("MS:1002049")) == null ? "null" : cvParam.getValue();
                     String specEValue =  (cvParam = cvParamMap.get("MS:1002052")) == null ? "null" : cvParam.getValue();
                     String eValue =  (cvParam = cvParamMap.get("MS:1002053")) == null ? "null" : cvParam.getValue();

                     Map<String, UserParam> userParamMap = getUserParamMap(sii.getUserParam());
                     UserParam userParam;
                     Integer isotopeError =  (userParam = userParamMap.get("IsotopeError")) == null ? null : Integer.parseInt(userParam.getValue());
                     
                     double precursorError = calculatedMassToCharge-(experimentalMassToCharge-Composition.ISOTOPE*isotopeError/charge);
                     
                     String fragmentationMethod = "N/A";
                     String protein = "N/A";
                     
                     out.println(specFile
                    		 +"\t"+specID
                    		 +"\t"+fragmentationMethod
                    		 +"\t"+calculatedMassToCharge
                    		 +"\t"+precursorError
                    		 +"\t"+charge
                    		 +"\t"+peptideSeq
                    		 +"\t"+protein
                    		 +"\t"+deNovoScore
                    		 +"\t"+rawScore
                    		 +"\t"+specEValue
                    		 +"\t"+eValue
                    		 );
                 } // end spectrum identification item
             } // end spectrum identification results
        }
        
        if(out != System.out)
        	out.close();
	}

	public static void main(String argv[]) throws Exception
	{
		MzMLAdapter.turnOffLogs();
		File mzidFile = new File(System.getProperty("user.home")+"/Research/Data/QCShew/test.mzid");
		MzIDParser parser = new MzIDParser(mzidFile);
		parser.writeToTSVFile(null);
	}	

}
