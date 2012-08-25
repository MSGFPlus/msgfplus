package edu.ucsd.msjava.mzid;

import java.io.File;
import java.util.List;

import edu.ucsd.msjava.mzml.MzMLAdapter;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisData;
import uk.ac.ebi.jmzidml.model.mzidml.DataCollection;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
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
	
	public static void main(String argv[]) throws Exception
	{
//		MzMLAdapter.turnOffLogs();
		File mzidFile = new File(System.getProperty("user.home")+"/Research/Data/QCShew/test.mzid");
		MzIDParser parser = new MzIDParser(mzidFile);
		parser.test();
	}	

}
