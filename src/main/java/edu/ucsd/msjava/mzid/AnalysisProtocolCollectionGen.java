package edu.ucsd.msjava.mzid;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.ucsd.msjava.msdbsearch.SearchParams;
import edu.ucsd.msjava.msutil.AminoAcidSet;

import uk.ac.ebi.jmzidml.model.mzidml.*;

public class AnalysisProtocolCollectionGen {
	private final SearchParams params;
	private final AminoAcidSet aaSet;
	private final Unimod unimod;
	private AnalysisProtocolCollection analysisProtocolCollection;
	private SpectrumIdentificationProtocol spectrumIdentificationProtocol;
	
	private Map<edu.ucsd.msjava.msutil.Modification, SearchModification> modMap;
	
	public AnalysisProtocolCollectionGen(SearchParams params, AminoAcidSet aaSet)
	{
		this.params = params;
		this.aaSet = aaSet;
		this.unimod = new Unimod();
		analysisProtocolCollection = new AnalysisProtocolCollection();
		generateSpectrumIdentificationProtocol();
	}

	public AnalysisProtocolCollection getAnalysisProtocolCollection()
	{
		return analysisProtocolCollection;
	}

	public SpectrumIdentificationProtocol getSpectrumIdentificationProtocol()
	{
		return spectrumIdentificationProtocol;
	}
	
	public SearchModification getSearchModification(edu.ucsd.msjava.msutil.Modification mod)
	{
		return modMap.get(mod);
	}
	
	private void generateSpectrumIdentificationProtocol()
	{
		spectrumIdentificationProtocol = new SpectrumIdentificationProtocol();
		analysisProtocolCollection.getSpectrumIdentificationProtocol().add(spectrumIdentificationProtocol);

        spectrumIdentificationProtocol.setId(Constants.siProtocolID);
        spectrumIdentificationProtocol.setAnalysisSoftware(Constants.msgfPlus);

        // SearchType
        Param searchTypeParam = new Param();
        searchTypeParam.setParam(Constants.makeCvParam("MS:1001083","ms-ms search"));
        spectrumIdentificationProtocol.setSearchType(searchTypeParam);

        // AdditionalSearchParams
        ParamList additionalSearchParams = new ParamList();
        List<CvParam> cvParamList = additionalSearchParams.getCvParam();
        cvParamList.add(Constants.makeCvParam("MS:1001211","parent mass type mono", Constants.psiCV));
        cvParamList.add(Constants.makeCvParam("MS:1001256","fragment mass type mono", Constants.psiCV));
        List<UserParam> userParamList = additionalSearchParams.getUserParam();
        userParamList.add(Constants.makeUserParam("MinIsotopeError", String.valueOf(params.getMin13C())));
        userParamList.add(Constants.makeUserParam("MaxIsotopeError", String.valueOf(params.getMax13C())));
        userParamList.add(Constants.makeUserParam("FragmentMethod", params.getActivationMethod().getName()));
        userParamList.add(Constants.makeUserParam("Instrument", params.getInstType().getName()));
        userParamList.add(Constants.makeUserParam("Protocol", params.getProtocol().getName()));
        userParamList.add(Constants.makeUserParam("NumTolerableTermini", String.valueOf(params.getNumTolerableTermini())));
        userParamList.add(Constants.makeUserParam("NumMatchesPerSpec", String.valueOf(params.getNumMatchesPerSpec())));
        // ModificationFile
        userParamList.add(Constants.makeUserParam("MinPepLength", String.valueOf(params.getMinPeptideLength())));
        userParamList.add(Constants.makeUserParam("MaxPepLength", String.valueOf(params.getMaxPeptideLength())));
        userParamList.add(Constants.makeUserParam("MinCharge", String.valueOf(params.getMinCharge())));
        userParamList.add(Constants.makeUserParam("MaxCharge", String.valueOf(params.getMaxCharge())));
        spectrumIdentificationProtocol.setAdditionalSearchParams(additionalSearchParams);
        
        ModificationParams modParams = getModificationParam();
        spectrumIdentificationProtocol.setModificationParams(modParams);

        // Enzymes
        Enzymes enzymes = new Enzymes();

        edu.ucsd.msjava.msutil.Enzyme enzyme = params.getEnzyme();
//        enzymes.setIndependent(false);
//        if(enzyme == null || enzyme == edu.ucsd.msjava.msutil.Enzyme.NOENZYME)
//        	enzymes.setIndependent(true);
//        else
//        	enzymes.setIndependent(false);

        if(enzyme != null)
        {
        	// Add enzyme
            List<Enzyme> enzymeList = enzymes.getEnzyme();
            Enzyme mzIdEnzyme = new Enzyme();
            mzIdEnzyme.setId(enzyme.getName());
            if(params.getNumTolerableTermini() == 0)
            	mzIdEnzyme.setSemiSpecific(false);
            else
            	mzIdEnzyme.setSemiSpecific(true);
            // Add name
            ParamList enzCvParams = new ParamList();
            String enzAcc = enzyme.getPSICvAccession();
            String enzName = enzyme.getDescription();
            
            if(enzAcc != null)	
                enzCvParams.getCvParam().add(Constants.makeCvParam(enzAcc, enzName, Constants.psiCV));
            else
            	enzCvParams.getUserParam().add(Constants.makeUserParam(enzName));
            mzIdEnzyme.setEnzymeName(enzCvParams);
            enzymeList.add(mzIdEnzyme);
        }

        spectrumIdentificationProtocol.setEnzymes(enzymes);

        // MassTable: skip
        
        // Fragment tolerance: N/A
        
        // Parent tolerance
        Tolerance parTol = new Tolerance();
        List<CvParam> parCvList = parTol.getCvParam();
        CvParam parCvPlus = Constants.getCvParamWithMassUnits(!params.getRightParentMassTolerance().isTolerancePPM());
        CvParam parCvMinus = Constants.getCvParamWithMassUnits(!params.getLeftParentMassTolerance().isTolerancePPM());
        parCvPlus.setAccession("MS:1001412");
        parCvPlus.setName("search tolerance plus value");
        parCvMinus.setAccession("MS:1001413");
        parCvMinus.setName("search tolerance minus value");
        parCvPlus.setValue(String.valueOf(params.getRightParentMassTolerance().getValue())); 
        parCvMinus.setValue(String.valueOf(params.getLeftParentMassTolerance().getValue()));
        parCvList.add(parCvPlus);
        parCvList.add(parCvMinus);
        spectrumIdentificationProtocol.setParentTolerance(parTol);

        // Threshold
        ParamList thrParamList = new ParamList();
        thrParamList.getCvParam().add(Constants.makeCvParam("MS:1001494","no threshold", Constants.psiCV));
        spectrumIdentificationProtocol.setThreshold(thrParamList);
	}
	
	public ModificationParams getModificationParam()
	{
		modMap = new HashMap<edu.ucsd.msjava.msutil.Modification, uk.ac.ebi.jmzidml.model.mzidml.SearchModification>();
		
        ModificationParams modParams = new ModificationParams();
        List<SearchModification> searchModList = modParams.getSearchModification();

        // fixed modifications
        for(edu.ucsd.msjava.msutil.Modification.Instance mod : aaSet.getModifications())
        {
    		String modName = mod.getModification().getName();
    		
            SearchModification searchMod = new SearchModification();

            searchMod.setFixedMod(true);
            searchMod.setMassDelta(mod.getModification().getMass());
    		
    		// set modification CV params
    		List<CvParam> modCvParamList = searchMod.getCvParam();
    		CvParam cvParam = new CvParam();
    		String unimodRecordID = unimod.getRecordID(modName);
    		if(unimodRecordID != null)	// exist in unimod
    		{
        		cvParam.setAccession(unimodRecordID);
        		cvParam.setCv(Constants.unimodCV);
        		cvParam.setName(modName);
    		}
    		else	// does not exist in Unimod
    		{
    			cvParam.setAccession("MS:1001460");	// unknown modification
    			cvParam.setName("unknown modification");
    		}
    		modCvParamList.add(cvParam);

    		// residue
            List<String> residueList = searchMod.getResidues();
            if(mod.getResidue() == '*')
            	residueList.add(".");
            else
            	residueList.add(String.valueOf(mod.getResidue()));

            // specificity rules
            SpecificityRules specificityRules = new SpecificityRules();
            List<CvParam> rules = specificityRules.getCvParam();
            if(mod.getLocation() != edu.ucsd.msjava.msutil.Modification.Location.Anywhere)
            	rules.add(Constants.makeCvParam("MS:1001055", "modification specificity rule", Constants.psiCV));
            if(mod.getLocation() == edu.ucsd.msjava.msutil.Modification.Location.N_Term || mod.getLocation() == edu.ucsd.msjava.msutil.Modification.Location.Protein_N_Term)
            	rules.add(Constants.makeCvParam("MS:1001189", "modification specificity N-term", Constants.psiCV));
            if(mod.getLocation() == edu.ucsd.msjava.msutil.Modification.Location.C_Term || mod.getLocation() == edu.ucsd.msjava.msutil.Modification.Location.Protein_C_Term)
            	rules.add(Constants.makeCvParam("MS:1001190", "modification specificity C-term", Constants.psiCV));
            searchMod.getSpecificityRules().add(specificityRules);
            
            searchModList.add(searchMod);
            
            modMap.put(mod.getModification(), searchMod);
        }

        return modParams;
	}
}