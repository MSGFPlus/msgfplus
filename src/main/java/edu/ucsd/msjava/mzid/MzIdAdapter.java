package edu.ucsd.msjava.mzid;

import java.util.List;
import java.util.Map;

import edu.ucsd.msjava.msgf.MSGFDBResultGenerator;
import edu.ucsd.msjava.params.ParamManager;


import uk.ac.ebi.jmzidml.model.mzidml.*;

public class MzIdAdapter {
	private SpectrumIdentificationProtocol searchParams;
	
	private List<DBSequence> dbSeqList;
	private List<Peptide> pepList;
	private List<PeptideEvidence> pepEviList;
	private SpectrumIdentificationList idList;
	
	public void addSearchParams(ParamManager params)
	{
		
	}
	
	private Map<String, Peptide> pepIDToPeptide;	// peptide id to peptide object
	private Map<String, DBSequence> protIDToProtein;	// db sequence id to db sequence object
	
	// set-up 
	public void addPeptideMatch(MSGFDBResultGenerator.DBMatch dbMatch)
	{
		
	}

	
	
}
