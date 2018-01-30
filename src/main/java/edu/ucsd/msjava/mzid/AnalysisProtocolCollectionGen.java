package edu.ucsd.msjava.mzid;

import edu.ucsd.msjava.msdbsearch.SearchParams;
import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Modification.Location;
import uk.ac.ebi.jmzidml.model.mzidml.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class AnalysisProtocolCollectionGen {
    private final SearchParams params;
    private final AminoAcidSet aaSet;
    private AnalysisProtocolCollection analysisProtocolCollection;
    private SpectrumIdentificationProtocol spectrumIdentificationProtocol;

    private Map<edu.ucsd.msjava.msutil.Modification, SearchModification> modMap;
    private Map<String, List<edu.ucsd.msjava.msutil.Modification>> fixedModMap;

    public AnalysisProtocolCollectionGen(SearchParams params, AminoAcidSet aaSet) {
        this.params = params;
        this.aaSet = aaSet;
        analysisProtocolCollection = new AnalysisProtocolCollection();
        modMap = new HashMap<edu.ucsd.msjava.msutil.Modification, uk.ac.ebi.jmzidml.model.mzidml.SearchModification>();
        fixedModMap = new HashMap<String, List<edu.ucsd.msjava.msutil.Modification>>();
        generateSpectrumIdentificationProtocol();
    }

    public AnalysisProtocolCollection getAnalysisProtocolCollection() {
        return analysisProtocolCollection;
    }

    public SpectrumIdentificationProtocol getSpectrumIdentificationProtocol() {
        return spectrumIdentificationProtocol;
    }

    public SearchModification getSearchModification(edu.ucsd.msjava.msutil.Modification mod) {
        return modMap.get(mod);
    }

    public List<edu.ucsd.msjava.msutil.Modification> getFixedModifications(char residue) {
        return fixedModMap.get(residue + "");
    }

    public List<edu.ucsd.msjava.msutil.Modification> getTerminalFixedModifications(char residue, Location location) {
        List<edu.ucsd.msjava.msutil.Modification> mods = new ArrayList<edu.ucsd.msjava.msutil.Modification>();
        if (location == Location.N_Term) {
            List<edu.ucsd.msjava.msutil.Modification> modList = fixedModMap.get("[" + residue);
            if (modList != null)
                mods.addAll(modList);

            modList = fixedModMap.get("[*");
            if (modList != null)
                mods.addAll(modList);
        } else if (location == Location.C_Term) {
            List<edu.ucsd.msjava.msutil.Modification> modList = fixedModMap.get(residue + "]");
            if (modList != null)
                mods.addAll(modList);

            modList = fixedModMap.get("*]");
            if (modList != null)
                mods.addAll(modList);
        }

        return mods;
    }

    private void generateSpectrumIdentificationProtocol() {
        spectrumIdentificationProtocol = new SpectrumIdentificationProtocol();
        analysisProtocolCollection.getSpectrumIdentificationProtocol().add(spectrumIdentificationProtocol);

        spectrumIdentificationProtocol.setId(Constants.siProtocolID);
        spectrumIdentificationProtocol.setAnalysisSoftware(Constants.msgfPlus);

        // SearchType
        Param searchTypeParam = new Param();
        searchTypeParam.setParam(Constants.makeCvParam("MS:1001083", "ms-ms search"));
        spectrumIdentificationProtocol.setSearchType(searchTypeParam);

        // Enzymes
        Enzymes enzymes = new Enzymes();
        edu.ucsd.msjava.msutil.Enzyme enzyme = params.getEnzyme();

        // AdditionalSearchParams
        ParamList additionalSearchParams = new ParamList();
        List<CvParam> cvParamList = additionalSearchParams.getCvParam();
        cvParamList.add(Constants.makeCvParam("MS:1001211", "parent mass type mono", Constants.psiCV));
        cvParamList.add(Constants.makeCvParam("MS:1001256", "fragment mass type mono", Constants.psiCV));
        List<UserParam> userParamList = additionalSearchParams.getUserParam();
        userParamList.add(Constants.makeUserParam("TargetDecoyApproach", String.valueOf(params.useTDA())));
        userParamList.add(Constants.makeUserParam("MinIsotopeError", String.valueOf(params.getMinIsotopeError())));
        userParamList.add(Constants.makeUserParam("MaxIsotopeError", String.valueOf(params.getMaxIsotopeError())));
        userParamList.add(Constants.makeUserParam("FragmentMethod", params.getActivationMethod().getName()));
        userParamList.add(Constants.makeUserParam("Instrument", params.getInstType().getName()));
        userParamList.add(Constants.makeUserParam("Protocol", params.getProtocol().getName()));
        int ntt = params.getNumTolerableTermini();
        if (enzyme == edu.ucsd.msjava.msutil.Enzyme.NoCleavage || enzyme == edu.ucsd.msjava.msutil.Enzyme.UnspecificCleavage)
            ntt = 0;
        userParamList.add(Constants.makeUserParam("NumTolerableTermini", String.valueOf(ntt)));
        userParamList.add(Constants.makeUserParam("NumMatchesPerSpec", String.valueOf(params.getNumMatchesPerSpec())));
        // ModificationFile
        userParamList.add(Constants.makeUserParam("MaxNumModifications", String.valueOf(aaSet.getMaxNumberOfVariableModificationsPerPeptide())));
        userParamList.add(Constants.makeUserParam("MinPepLength", String.valueOf(params.getMinPeptideLength())));
        userParamList.add(Constants.makeUserParam("MaxPepLength", String.valueOf(params.getMaxPeptideLength())));
        userParamList.add(Constants.makeUserParam("MinCharge", String.valueOf(params.getMinCharge())));
        userParamList.add(Constants.makeUserParam("MaxCharge", String.valueOf(params.getMaxCharge())));
        userParamList.add(Constants.makeUserParam("ChargeCarrierMass", String.valueOf(params.getChargeCarrierMass())));
        spectrumIdentificationProtocol.setAdditionalSearchParams(additionalSearchParams);

        if (!aaSet.getModifications().isEmpty()) {
            ModificationParams modParams = getModificationParam();
            spectrumIdentificationProtocol.setModificationParams(modParams);
        }

//        enzymes.setIndependent(false);
//        if(enzyme == null || enzyme == edu.ucsd.msjava.msutil.Enzyme.NOENZYME)
//        	enzymes.setIndependent(true);
//        else
//        	enzymes.setIndependent(false);

        if (enzyme != null) {
            // Add enzyme
            List<Enzyme> enzymeList = enzymes.getEnzyme();
            Enzyme mzIdEnzyme = new Enzyme();
            mzIdEnzyme.setId(enzyme.getName());
            if (ntt == 2)
                mzIdEnzyme.setSemiSpecific(false);
            else
                mzIdEnzyme.setSemiSpecific(true);
            mzIdEnzyme.setMissedCleavages(1000);
            // Add name
            ParamList enzCvParams = new ParamList();
            String enzAcc = enzyme.getPSICvAccession();
            String enzName = enzyme.getDescription();

            if (enzAcc != null)
                enzCvParams.getCvParam().add(Constants.makeCvParam(enzAcc, enzName, Constants.psiCV));
            else
                enzCvParams.getUserParam().add(Constants.makeUserParam(enzName));
            mzIdEnzyme.setEnzymeName(enzCvParams);
            enzymeList.add(mzIdEnzyme);
        }

        spectrumIdentificationProtocol.setEnzymes(enzymes);

        // MassTable: Only output custom residues
        ArrayList<Character> defaultResidues = AminoAcidSet.getStandardAminoAcidSet().getResidueListWithoutMods();
        ArrayList<Character> usedResidues = this.aaSet.getResidueListWithoutMods();
        if (usedResidues.size() > defaultResidues.size()) {
            MassTable massTable = new MassTable();
            massTable.setId(Constants.massTableID);
            massTable.setName("Custom Residues");
            massTable.getMsLevel().add(1);
            massTable.getMsLevel().add(2);
            for (Character c : usedResidues) {
                // Only outputting the masses of non-standard residues
                if (!defaultResidues.contains(c)) {
                    AminoAcid aa = this.aaSet.getAminoAcid(c);
                    Residue residue = new Residue();
                    residue.setCode(aa.getResidueStr());
                    residue.setMass(aa.getMass());
                    massTable.getResidue().add(residue);
                }
            }
            spectrumIdentificationProtocol.getMassTable().add(massTable);
        }

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
        thrParamList.getCvParam().add(Constants.makeCvParam("MS:1001494", "no threshold", Constants.psiCV));
        spectrumIdentificationProtocol.setThreshold(thrParamList);
    }

    public ModificationParams getModificationParam() {
        ModificationParams modParams = new ModificationParams();
        List<SearchModification> searchModList = modParams.getSearchModification();

        // fixed modifications
        for (edu.ucsd.msjava.msutil.Modification.Instance mod : aaSet.getModifications()) {
            String modName = mod.getModification().getName();

            SearchModification searchMod = new SearchModification();

            searchMod.setFixedMod(mod.isFixedModification());
            searchMod.setMassDelta(mod.getModification().getMass());

            // set modification CV params
            List<CvParam> modCvParamList = searchMod.getCvParam();
            CvParam cvParam = new CvParam();
            String unimodRecordID = Unimod.getUnimod().getRecordID(modName);
            if (unimodRecordID != null)    // exist in unimod
            {
                cvParam.setAccession(unimodRecordID);
                cvParam.setCv(Constants.unimodCV);
                cvParam.setName(modName);
            } else    // does not exist in Unimod
            {
                cvParam.setAccession("MS:1001460");    // unknown modification
                cvParam.setCv(Constants.psiCV);
                cvParam.setName("unknown modification");
                cvParam.setValue(modName);
            }
            modCvParamList.add(cvParam);

            // residue
            List<String> residueList = searchMod.getResidues();
            if (mod.getResidue() == '*')
                residueList.add(".");
            else
                residueList.add(String.valueOf(mod.getResidue()));

            if (mod.getLocation() != edu.ucsd.msjava.msutil.Modification.Location.Anywhere) {
                // specificity rules
                SpecificityRules specificityRules = new SpecificityRules();
                List<CvParam> rules = specificityRules.getCvParam();
                if (mod.getLocation() == edu.ucsd.msjava.msutil.Modification.Location.N_Term)
                    rules.add(Constants.makeCvParam("MS:1001189", "modification specificity N-term", Constants.psiCV));
                else if (mod.getLocation() == edu.ucsd.msjava.msutil.Modification.Location.Protein_N_Term)
                    rules.add(Constants.makeCvParam("MS:1002057", "modification specificity protein N-term", Constants.psiCV));
                else if (mod.getLocation() == edu.ucsd.msjava.msutil.Modification.Location.C_Term)
                    rules.add(Constants.makeCvParam("MS:1001190", "modification specificity C-term", Constants.psiCV));
                else if (mod.getLocation() == edu.ucsd.msjava.msutil.Modification.Location.Protein_C_Term)
                    rules.add(Constants.makeCvParam("MS:1002058", "modification specificity protein C-term", Constants.psiCV));
                searchMod.getSpecificityRules().add(specificityRules);
            }

            searchModList.add(searchMod);

            modMap.put(mod.getModification(), searchMod);
            if (mod.isFixedModification()) {
                String modKey = getModKey(mod.getResidue(), mod.getLocation());
                List<edu.ucsd.msjava.msutil.Modification> fixedMods = fixedModMap.get(modKey);
                if (fixedMods == null) {
                    fixedMods = new ArrayList<edu.ucsd.msjava.msutil.Modification>();
                    fixedModMap.put(modKey, fixedMods);
                }
                fixedMods.add(mod.getModification());
            }
        }

        return modParams;
    }

    private String getModKey(char residue, Location location) {
        StringBuffer modKey = new StringBuffer();

        if (location == Location.N_Term)
            modKey.append('[');

        modKey.append(residue);

        if (location == Location.C_Term)
            modKey.append(']');

        return modKey.toString();
    }
}