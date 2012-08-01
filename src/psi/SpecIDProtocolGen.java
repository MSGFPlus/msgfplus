//package psi;
//
//import java.util.List;
//
//import msutil.AminoAcidSet;
//
//import uk.ac.ebi.jmzidml.model.mzidml.*;
//
//public class SpecIDProtocolGen {
//	private msutil.AminoAcidSet aaSet;
//	private AnalysisSoftware analysisSoftware;
//	private msutil.Enzyme enzyme;
//	private int ntt;
//	private msgf.Tolerance leftPrecursorMassTolerance;
//	private msgf.Tolerance rightPrecursorMassTolerance;
//	private int fragmentMethodID;
//	private int instrumentID;
//	private int numAllowedC13;
//	
//	public SpecIDProtocolGen(
//			AminoAcidSet aaSet,
//			AnalysisSoftware analysisSoftware,
//			msutil.Enzyme enzyme,
//			int ntt,					// number of tolerable termini
//			msgf.Tolerance leftPrecursorMassTolerance,
//			msgf.Tolerance rightPrecursorMassTolerance,
//			int numAllowedC13
//			)
//	{
//		this.aaSet = aaSet;
//		this.analysisSoftware = analysisSoftware;
//		this.enzyme = enzyme;
//		this.ntt = ntt;
//		this.leftPrecursorMassTolerance = leftPrecursorMassTolerance;
//		this.rightPrecursorMassTolerance = rightPrecursorMassTolerance;
//		this.numAllowedC13 = numAllowedC13;
//	}
//	
//	public SpectrumIdentificationProtocol genSpectrumIdentificationProtocol()
//	{
//		SpectrumIdentificationProtocol siProtocol = new SpectrumIdentificationProtocol();
//        siProtocol.setId(Constants.siProtocolID);
//        siProtocol.setAnalysisSoftware(analysisSoftware);
//
//        // SearchType
//        Param searchTypeParam = new Param();
//        searchTypeParam.setParam(Constants.makeCvParam("MS:1001083","ms-ms search", Constants.psiCV));
//        siProtocol.setSearchType(searchTypeParam);
//
//        // AdditionalSearchParams
//        // TODO: add fragmentation, instrument, protocol, numAllowedC13
//        
//        ParamList additionalSearchParams = new ParamList();
//        List<CvParam> cvParamList = additionalSearchParams.getCvParam();
//        cvParamList.add(Constants.makeCvParam("MS:1001211","parent mass type mono", Constants.psiCV));
//        cvParamList.add(Constants.makeCvParam("MS:1001256","fragment mass type mono", Constants.psiCV));
//        siProtocol.setAdditionalSearchParams(additionalSearchParams);
//
//        List<UserParam> usrParamList = additionalSearchParams.getUserParam();
//        usrParamList.add(Constants.makeUserParam(name));
//        
//        // ModificationParams
//        ModificationParams modParams = getModificationParam();
//        siProtocol.setModificationParams(modParams);
//
//        // Enzymes
//        Enzymes enzymes = new Enzymes();
//
//        if(enzyme == null)
//        	enzymes.setIndependent(true);
//        else
//        {
//        	enzymes.setIndependent(false);
//        	
//        	// Add enzyme
//            List<Enzyme> enzymeList = enzymes.getEnzyme();
//            Enzyme mzIdEnzyme = new Enzyme();
//            if(ntt == 0)
//            	mzIdEnzyme.setSemiSpecific(false);
//            else
//            	mzIdEnzyme.setSemiSpecific(true);
//            // Add name
//            ParamList enzCvParams = new ParamList();
//            String enzAcc = enzyme.getPSICvAccession();
//            String enzName = enzyme.getName();
//            
//            if(enzAcc != null)	
//                enzCvParams.getCvParam().add(Constants.makeCvParam(enzAcc, enzName, Constants.psiCV));
//            else
//            	enzCvParams.getUserParam().add(Constants.makeUserParam(enzName));
//            mzIdEnzyme.setEnzymeName(enzCvParams);
//            enzymeList.add(mzIdEnzyme);
//        }
//
//        siProtocol.setEnzymes(enzymes);
//
//        // MassTable: skip
//        
//        // Fragment tolerance: no
//        
//        // Parent tolerance
//
//        Tolerance parTol = new Tolerance();
//
//        List<CvParam> parCvList = parTol.getCvParam();
//        CvParam parCvPlus = Constants.getCvParamWithMassUnits(!this.rightPrecursorMassTolerance.isTolerancePPM());
//        CvParam parCvMinus = Constants.getCvParamWithMassUnits(!this.leftPrecursorMassTolerance.isTolerancePPM());
//
//        parCvPlus.setAccession("MS:1001412");
//        parCvPlus.setName("search tolerance plus value");
//        parCvMinus.setAccession("MS:1001413");
//        parCvMinus.setName("search tolerance minus value");
//        parCvPlus.setValue(String.valueOf(this.rightPrecursorMassTolerance.getValue())); 
//        parCvMinus.setValue(String.valueOf(this.leftPrecursorMassTolerance.getValue()));
//        parCvList.add(parCvPlus);
//        parCvList.add(parCvMinus);
//
//        siProtocol.setParentTolerance(parTol);
//
//        // Threshold
//        ParamList thrParamList = new ParamList();
//        thrParamList.getCvParam().add(Constants.makeCvParam("MS:1001494","no threshold", Constants.psiCV));
//        //<cvParam accession="MS:1001494" name="no threshold" cvRef="PSI-MS" />
//        siProtocol.setThreshold(thrParamList);
//        
//        return siProtocol;
//	}
//	
//	public ModificationParams getModificationParam()
//	{
//        ModificationParams modParams = new ModificationParams();
//        List<SearchModification> searchModList = modParams.getSearchModification();
//
//        // fixed modifications
//        for(msutil.Modification.Instance mod : aaSet.getModifications())
//        {
//    		String modName = mod.getModification().getName();
//    		
//            SearchModification searchMod = new SearchModification();
//
//            searchMod.setFixedMod(true);
//            searchMod.setMassDelta(mod.getModification().getMass());
//    		
//    		// set modification CV params
//    		List<CvParam> modCvParamList = searchMod.getCvParam();
//    		CvParam cvParam = new CvParam();
//    		int unimodRecordID = Constants.unimod.getRecordID(modName);
//    		if(unimodRecordID > 0)	// exist in unimod
//    		{
//        		cvParam.setAccession("UNIMOD:" + unimodRecordID);
//        		cvParam.setCv(Constants.unimodCV);
//        		cvParam.setName(modName);
//    		}
//    		else	// does not exist in Unimod
//    		{
//    			cvParam.setAccession("MS:1001460");	// unknown modification
//    			cvParam.setName("unknown modification");
//    		}
//    		modCvParamList.add(cvParam);
//
//    		// residue
//            List<String> residueList = searchMod.getResidues();
//            residueList.add(String.valueOf(mod.getResidue()));	// TODO: what if residue='*'
//
//            // specificity rules
//            SpecificityRules specificityRules = new SpecificityRules();
//            List<CvParam> rules = specificityRules.getCvParam();
//            rules.add(Constants.makeCvParam("MS:1001055", "modification specificity rule", Constants.psiCV));
//            if(mod.getLocation() == msutil.Modification.Location.N_Term || mod.getLocation() == msutil.Modification.Location.Protein_N_Term)
//            	rules.add(Constants.makeCvParam("MS:1001189", "modification specificity N-term", Constants.psiCV));
//            if(mod.getLocation() == msutil.Modification.Location.C_Term || mod.getLocation() == msutil.Modification.Location.Protein_C_Term)
//            	rules.add(Constants.makeCvParam("MS:1001190", "modification specificity C-term", Constants.psiCV));
//            searchModList.add(searchMod);
//        }
//
//        return modParams;
//	}
//}
